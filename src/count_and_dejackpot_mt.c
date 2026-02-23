/*
 * count_and_dejackpot_mt.c  –  optimised multi-threaded C port
 *
 * Identical logic to count_and_dejackpot.c, with POSIX thread parallelism.
 * Use --threads N to control worker count (default: 1).
 *
 * Read structure:
 *   [outer5] [umi5] [inner5] [barcode] [inner3] [umi3] [outer3]
 *
 * ── COMPILATION ──────────────────────────────────────────────────────────────
 *
 * Dependencies:
 *   1. zlib       – gzip FASTQ support
 *                   Ubuntu/Debian : sudo apt install zlib1g-dev
 *                   macOS (brew)  : brew install zlib
 *
 *   2. uthash.h   – single-header hash table library (no build step)
 *                   Homepage: https://troydhanson.github.io/uthash/
 *                   Download one file:
 *                     wget https://raw.githubusercontent.com/troydhanson/uthash/master/src/uthash.h
 *                   Place uthash.h next to this file (or anywhere on your include path).
 *
 *   3. pthreads   – POSIX threads (included with glibc; link with -lpthread)
 *
 * Recommended build commands:
 *   gcc  -O2 -o count_and_dejackpot_mt count_and_dejackpot_mt.c -lz -lm -lpthread
 *   clang -O2 -o count_and_dejackpot_mt count_and_dejackpot_mt.c -lz -lm -lpthread
 *
 * ── THREADING DESIGN ─────────────────────────────────────────────────────────
 *
 * Producer / consumer with batched queue:
 *
 *   Main thread (producer)
 *     Reads FASTQ records and accumulates them into batches of BATCH_SIZE.
 *     Each full batch is pushed into the shared queue as one unit, amortising
 *     mutex overhead across BATCH_SIZE records.
 *
 *   Worker threads (consumers)
 *     Each worker pops one batch at a time from the queue, processes every
 *     record in it, and accumulates parsed results into its own thread-local
 *     TupleEntry hash table.  Failed reads are written to a per-thread memory
 *     buffer and flushed to the shared fail file in bulk, eliminating per-read
 *     mutex contention on the fail path.
 *
 *   Merge step (serial, after all workers finish)
 *     The main thread folds all per-worker TupleEntry tables into one global
 *     table.  No locking needed at this point.
 *
 *   Everything after the merge is identical to the single-threaded version.
 *
 * ── PERFORMANCE NOTES ────────────────────────────────────────────────────────
 *
 * Key optimisations vs the naive multi-threaded version:
 *   1. Batched queue   – BATCH_SIZE records per slot → mutex ops reduced by
 *                        BATCH_SIZE×.  Dominant when n_threads is high.
 *   2. Live-byte copy  – queue push/pop copy only strlen(field)+1 bytes, not
 *                        the full MAX_LINE buffer.  Cuts per-record memory
 *                        traffic by ~60× for typical FLASH-merged reads.
 *   3. Fail buffer     – per-thread; flushes to disk only when full or at exit.
 *                        Removes the fail-file mutex from the hot path entirely.
 *   4. Direct key build – build_tuple_key() uses memcpy with pre-computed
 *                        lengths from ParseResult instead of snprintf.
 *   5. 1 MB gzFile buffer – fewer syscalls for both plain and gzipped input.
 *
 * For gzipped input, gzip decompression is single-threaded (the producer)
 * and often limits throughput regardless of worker count.  For plain FASTQ
 * on fast storage (NVMe / ramdisk) or piped input, worker threads can provide
 * close to linear speedup up to the number of physical cores.
 *
 * ── FUZZY MATCHING NOTE ──────────────────────────────────────────────────────
 *
 * Substitutions only (Hamming distance); see count_and_dejackpot.c for rationale.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <pthread.h>
#include <zlib.h>
#include "uthash.h"

/* ── constants ──────────────────────────────────────────────────────────── */

#define MAX_LINE    16384  /* max bytes in a single FASTQ line (incl. newline) */
#define MAX_UMI       256  /* max UMI sequence length                          */
#define MAX_BC        256  /* max barcode sequence length                      */

#define MAX_TUPLE_KEY  (MAX_UMI + 1 + MAX_BC + 1 + MAX_UMI + 1)
#define MAX_UMI_KEY    (MAX_UMI + 1 + MAX_UMI + 1)

/* Queue / threading tuning – safe to adjust at compile time */
#define BATCH_SIZE          32    /* FASTQ records per queue slot             */
#define MIN_QUEUE_BATCHES    4    /* minimum queue depth (in batches)         */
#define MAX_QUEUE_BATCHES   64    /* maximum queue depth (in batches)         */
#define FAIL_BUF_BYTES   (256 * 1024)  /* per-thread fail-record buffer       */

/* ── shared data structures ─────────────────────────────────────────────── */

typedef struct {
    const char *inner5;  int inner5_len;
    const char *inner3;  int inner3_len;
    const char *outer5;  int outer5_len;
    const char *outer3;  int outer3_len;
    int barcode_length;
    int min_umi, max_umi;
    int max_errors;
} ParseArgs;

/* Lengths stored alongside strings to avoid strlen in the hot path */
typedef struct {
    char umi5   [MAX_UMI]; int u5len;
    char barcode[MAX_BC];  int bclen;
    char umi3   [MAX_UMI]; int u3len;
} ParseResult;

typedef struct TupleEntry {
    char           key[MAX_TUPLE_KEY];
    long           count;
    UT_hash_handle hh;
} TupleEntry;

typedef struct UmiPairEntry {
    char           key[MAX_UMI_KEY];
    UT_hash_handle hh;
} UmiPairEntry;

typedef struct BarcodeSummary {
    char           barcode[MAX_BC];
    long           total_reads;
    long           num_unique_umis;
    UmiPairEntry  *umi_set;
    UT_hash_handle hh;
} BarcodeSummary;

/* ── FASTQ record and batch ─────────────────────────────────────────────── */

/* Stores live string lengths so push/pop can copy only the relevant bytes. */
typedef struct {
    char header[MAX_LINE]; int hlen;
    char seq   [MAX_LINE]; int slen;
    char qual  [MAX_LINE]; int qlen;
} FastqRecord;

typedef struct {
    FastqRecord recs[BATCH_SIZE];
    int         count;
} RecordBatch;

/* ── work queue (batched circular buffer) ───────────────────────────────── */

typedef struct {
    RecordBatch    *slots;
    int             capacity;   /* depth in batches                           */
    int             head, tail, count;
    int             done;
    pthread_mutex_t mutex;
    pthread_cond_t  not_empty;
    pthread_cond_t  not_full;
} WorkQueue;

static void queue_init(WorkQueue *q, int capacity)
{
    q->slots = malloc((size_t)capacity * sizeof(RecordBatch));
    if (!q->slots) { perror("malloc"); exit(1); }
    q->capacity = capacity;
    q->head = q->tail = q->count = q->done = 0;
    pthread_mutex_init(&q->mutex,     NULL);
    pthread_cond_init (&q->not_empty, NULL);
    pthread_cond_init (&q->not_full,  NULL);
}

static void queue_destroy(WorkQueue *q)
{
    free(q->slots);
    pthread_mutex_destroy(&q->mutex);
    pthread_cond_destroy (&q->not_empty);
    pthread_cond_destroy (&q->not_full);
}

/*
 * Copy only live bytes (hlen/slen/qlen + 1) from src to dst.
 * For a ~400 bp read this copies ~850 bytes instead of 49 KB.
 */
static inline void copy_record(FastqRecord *dst, const FastqRecord *src)
{
    dst->hlen = src->hlen; memcpy(dst->header, src->header, (size_t)src->hlen + 1);
    dst->slen = src->slen; memcpy(dst->seq,    src->seq,    (size_t)src->slen + 1);
    dst->qlen = src->qlen; memcpy(dst->qual,   src->qual,   (size_t)src->qlen + 1);
}

/* Push one batch into the queue (blocks if full). */
static void queue_push_batch(WorkQueue *q, const RecordBatch *batch)
{
    pthread_mutex_lock(&q->mutex);
    while (q->count == q->capacity)
        pthread_cond_wait(&q->not_full, &q->mutex);

    RecordBatch *slot = &q->slots[q->tail];
    slot->count = batch->count;
    for (int i = 0; i < batch->count; i++)
        copy_record(&slot->recs[i], &batch->recs[i]);

    q->tail = (q->tail + 1) % q->capacity;
    q->count++;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mutex);
}

/*
 * Pop one batch from the queue into `out` (blocks until available or done).
 * Returns 1 on success, 0 when queue is drained and done.
 */
static int queue_pop_batch(WorkQueue *q, RecordBatch *out)
{
    pthread_mutex_lock(&q->mutex);
    while (q->count == 0 && !q->done)
        pthread_cond_wait(&q->not_empty, &q->mutex);

    if (q->count == 0) {
        pthread_mutex_unlock(&q->mutex);
        return 0;
    }

    const RecordBatch *slot = &q->slots[q->head];
    out->count = slot->count;
    for (int i = 0; i < slot->count; i++)
        copy_record(&out->recs[i], &slot->recs[i]);

    q->head = (q->head + 1) % q->capacity;
    q->count--;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mutex);
    return 1;
}

static void queue_finish(WorkQueue *q)
{
    pthread_mutex_lock(&q->mutex);
    q->done = 1;
    pthread_cond_broadcast(&q->not_empty);
    pthread_mutex_unlock(&q->mutex);
}

/* ── per-thread fail-record buffer ─────────────────────────────────────── */

typedef struct {
    char *data;
    int   used;
} FailBuffer;

static void fail_buf_flush(FailBuffer *fb, FILE *fp, pthread_mutex_t *mu)
{
    if (fb->used == 0) return;
    pthread_mutex_lock(mu);
    fwrite(fb->data, 1, (size_t)fb->used, fp);
    pthread_mutex_unlock(mu);
    fb->used = 0;
}

/* Append one failed record; flushes to disk when the buffer is nearly full. */
static void fail_buf_write(FailBuffer *fb,
                            const char *hdr, int hlen,
                            const char *seq, int slen,
                            const char *qual, int qlen,
                            FILE *fp, pthread_mutex_t *mu)
{
    /* hdr\n + seq\n + +\n + qual\n = hlen+1 + slen+1 + 2 + qlen+1 = +5 */
    int need = hlen + slen + qlen + 5;
    if (fb->used + need > FAIL_BUF_BYTES)
        fail_buf_flush(fb, fp, mu);
    /* If one record alone exceeds the buffer, write directly */
    if (need > FAIL_BUF_BYTES) {
        pthread_mutex_lock(mu);
        fprintf(fp, "%s\n%s\n+\n%s\n", hdr, seq, qual);
        pthread_mutex_unlock(mu);
        return;
    }
    char *p = fb->data + fb->used;
    memcpy(p, hdr,  (size_t)hlen); p += hlen; *p++ = '\n';
    memcpy(p, seq,  (size_t)slen); p += slen; *p++ = '\n';
    *p++ = '+'; *p++ = '\n';
    memcpy(p, qual, (size_t)qlen); p += qlen; *p++ = '\n';
    fb->used = (int)(p - fb->data);
}

/* ── worker thread state ────────────────────────────────────────────────── */

typedef struct {
    WorkQueue       *queue;
    const ParseArgs *pa;
    TupleEntry      *local_tuples;
    long             n_parsed;
    long             n_failed;
    FILE            *fail_fp;
    pthread_mutex_t *fail_mutex;
} WorkerArgs;

/* ── fuzzy matching (stateless, thread-safe) ────────────────────────────── */

static inline int count_subst(const char *text, int pos,
                               const char *pattern, int plen, int limit)
{
    int err = 0;
    for (int j = 0; j < plen; j++)
        if (text[pos + j] != pattern[j])
            if (++err > limit) return err;
    return err;
}

static int find_core(const char *seq, int slen, const ParseArgs *pa,
                     int *i5_start, int *i5_end,
                     int *i3_start, int *i3_end)
{
    int span = pa->inner5_len + pa->barcode_length + pa->inner3_len;
    if (span > slen) return -1;

    int best_pos = -1;
    int best_err = 2 * pa->max_errors + 1;

    for (int pos = 0; pos <= slen - span; pos++) {
        int e5 = count_subst(seq, pos, pa->inner5, pa->inner5_len, pa->max_errors);
        if (e5 > pa->max_errors) continue;

        int i3_pos = pos + pa->inner5_len + pa->barcode_length;
        int e3 = count_subst(seq, i3_pos, pa->inner3, pa->inner3_len, pa->max_errors);
        if (e3 > pa->max_errors) continue;

        int total = e5 + e3;
        if (total < best_err) {
            best_err = total;
            best_pos = pos;
            if (total == 0) break;
        }
    }

    if (best_pos < 0) return -1;

    *i5_start = best_pos;
    *i5_end   = best_pos + pa->inner5_len;
    *i3_start = *i5_end  + pa->barcode_length;
    *i3_end   = *i3_start + pa->inner3_len;
    return 0;
}

static int fuzzy_find(const char *text, int tlen,
                      const char *pattern, int plen,
                      int max_errors, int *match_end)
{
    if (plen <= 0 || plen > tlen) return -1;
    int best_pos = -1, best_err = max_errors + 1;

    for (int i = 0; i <= tlen - plen; i++) {
        int err = count_subst(text, i, pattern, plen, max_errors);
        if (err < best_err) {
            best_err = err;
            best_pos = i;
            if (err == 0) break;
        }
    }

    if (best_pos < 0 || best_err > max_errors) return -1;
    if (match_end) *match_end = best_pos + plen;
    return best_pos;
}

/* ── read parsing (stateless, thread-safe) ──────────────────────────────── */

static int parse_read(const char *seq, const ParseArgs *pa, ParseResult *out)
{
    int slen = (int)strlen(seq);
    int i5_start, i5_end, i3_start, i3_end;

    if (find_core(seq, slen, pa, &i5_start, &i5_end, &i3_start, &i3_end) < 0)
        return -1;

    int u5_start, u5_end = i5_start;
    if (pa->outer5) {
        int o5_end;
        if (fuzzy_find(seq, i5_start, pa->outer5, pa->outer5_len,
                       pa->max_errors, &o5_end) < 0)
            return -1;
        u5_start = o5_end;
    } else {
        u5_start = 0;
    }
    int u5_len = u5_end - u5_start;
    if (u5_len < pa->min_umi || u5_len > pa->max_umi) return -1;

    int u3_start = i3_end, u3_end;
    if (pa->outer3) {
        int o3_start = fuzzy_find(seq + i3_end, slen - i3_end,
                                  pa->outer3, pa->outer3_len, pa->max_errors, NULL);
        if (o3_start < 0) return -1;
        u3_end = i3_end + o3_start;
    } else {
        u3_end = slen;
    }
    int u3_len = u3_end - u3_start;
    if (u3_len < pa->min_umi || u3_len > pa->max_umi) return -1;

    if (u5_len >= MAX_UMI || pa->barcode_length >= MAX_BC || u3_len >= MAX_UMI)
        return -1;

    memcpy(out->umi5,    seq + u5_start,      (size_t)u5_len);            out->umi5[u5_len] = '\0';
    memcpy(out->barcode, seq + i5_end, (size_t)pa->barcode_length); out->barcode[pa->barcode_length] = '\0';
    memcpy(out->umi3,    seq + u3_start,      (size_t)u3_len);            out->umi3[u3_len] = '\0';

    /* Store lengths to avoid strlen in caller */
    out->u5len = u5_len;
    out->bclen = pa->barcode_length;
    out->u3len = u3_len;
    return 0;
}

/* ── fast tuple key construction ────────────────────────────────────────── */

/*
 * Build "umi5\tbarcode\tumi3" into key[0..keymax) using memcpy with
 * pre-computed lengths.  Returns total key length (excl. '\0'), or -1 if
 * the key would not fit.  Avoids the formatting overhead of snprintf.
 */
static inline int build_tuple_key(char *key, int keymax,
                                   const char *u5, int u5l,
                                   const char *bc, int bcl,
                                   const char *u3, int u3l)
{
    int total = u5l + 1 + bcl + 1 + u3l;
    if (total >= keymax) return -1;
    char *p = key;
    memcpy(p, u5, (size_t)u5l); p += u5l; *p++ = '\t';
    memcpy(p, bc, (size_t)bcl); p += bcl; *p++ = '\t';
    memcpy(p, u3, (size_t)u3l); p += u3l; *p   = '\0';
    return total;
}

/* ── worker thread ──────────────────────────────────────────────────────── */

static void *worker_thread(void *arg)
{
    WorkerArgs *wa = arg;

    /* Heap-allocate large working buffers to keep the stack frame small */
    RecordBatch *batch = malloc(sizeof(RecordBatch));
    FailBuffer   fb    = { malloc(FAIL_BUF_BYTES), 0 };
    if (!batch || !fb.data) { perror("malloc"); exit(1); }

    while (queue_pop_batch(wa->queue, batch)) {
        for (int i = 0; i < batch->count; i++) {
            const FastqRecord *r = &batch->recs[i];
            ParseResult pr;
            int parsed = 0;

            if (parse_read(r->seq, wa->pa, &pr) == 0) {
                char key[MAX_TUPLE_KEY];
                int klen = build_tuple_key(key, (int)sizeof(key),
                                           pr.umi5,    pr.u5len,
                                           pr.barcode, pr.bclen,
                                           pr.umi3,    pr.u3len);
                if (klen > 0) {
                    TupleEntry *e = NULL;
                    HASH_FIND_STR(wa->local_tuples, key, e);
                    if (!e) {
                        e = malloc(sizeof(TupleEntry));
                        if (!e) { perror("malloc"); exit(1); }
                        memcpy(e->key, key, (size_t)klen + 1);
                        e->count = 0;
                        HASH_ADD_STR(wa->local_tuples, key, e);
                    }
                    e->count++;
                    wa->n_parsed++;
                    parsed = 1;
                }
            }

            if (!parsed) {
                wa->n_failed++;
                fail_buf_write(&fb,
                               r->header, r->hlen,
                               r->seq,    r->slen,
                               r->qual,   r->qlen,
                               wa->fail_fp, wa->fail_mutex);
            }
        }
    }

    fail_buf_flush(&fb, wa->fail_fp, wa->fail_mutex);
    free(fb.data);
    free(batch);
    return NULL;
}

/* ── merge per-worker tables into global ────────────────────────────────── */

static void merge_tables(TupleEntry **global, TupleEntry *local)
{
    TupleEntry *e, *tmp;
    HASH_ITER(hh, local, e, tmp) {
        HASH_DEL(local, e);
        TupleEntry *g = NULL;
        HASH_FIND_STR(*global, e->key, g);
        if (g) {
            g->count += e->count;
            free(e);
        } else {
            HASH_ADD_STR(*global, key, e);
        }
    }
}

/* ── Gini coefficient ───────────────────────────────────────────────────── */

static int cmp_long_asc(const void *a, const void *b)
{
    long x = *(const long *)a, y = *(const long *)b;
    return (x > y) - (x < y);
}

static double gini(const long *values, int n)
{
    if (n <= 0) return NAN;
    long *s = malloc((size_t)n * sizeof(long));
    if (!s) { perror("malloc"); exit(1); }
    memcpy(s, values, (size_t)n * sizeof(long));
    qsort(s, (size_t)n, sizeof(long), cmp_long_asc);

    long total = 0;
    for (int i = 0; i < n; i++) total += s[i];
    if (total == 0) { free(s); return 0.0; }

    double num = 0.0;
    for (int i = 0; i < n; i++)
        num += (double)s[i] * (2.0 * (i + 1) - n - 1);
    free(s);
    return num / ((double)n * (double)total);
}

/* ── FASTQ I/O ──────────────────────────────────────────────────────────── */

/*
 * Fills r->header/seq/qual and pre-computes r->hlen/slen/qlen so callers
 * never need strlen on the hot path.
 */
static int read_fastq_record(gzFile fp, FastqRecord *r)
{
    char plus[MAX_LINE];
    if (!gzgets(fp, r->header, MAX_LINE) || r->header[0] == '\0') return 0;
    if (!gzgets(fp, r->seq,    MAX_LINE)) return 0;
    if (!gzgets(fp, plus,      MAX_LINE)) return 0;
    if (!gzgets(fp, r->qual,   MAX_LINE)) return 0;
    r->hlen = (int)strcspn(r->header, "\r\n"); r->header[r->hlen] = '\0';
    r->slen = (int)strcspn(r->seq,    "\r\n"); r->seq   [r->slen] = '\0';
    r->qlen = (int)strcspn(r->qual,   "\r\n"); r->qual  [r->qlen] = '\0';
    return 1;
}

/* ── sort comparator ────────────────────────────────────────────────────── */

static int cmp_bc_desc(const void *a, const void *b)
{
    const BarcodeSummary *x = *(const BarcodeSummary * const *)a;
    const BarcodeSummary *y = *(const BarcodeSummary * const *)b;
    if (y->num_unique_umis > x->num_unique_umis) return  1;
    if (y->num_unique_umis < x->num_unique_umis) return -1;
    return 0;
}

/* ── usage ──────────────────────────────────────────────────────────────── */

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s --in <fastq[.gz]> --out <counts.tsv> --fail <unparsed.fastq>\n"
        "          --inner5 <seq> --inner3 <seq> --barcode-length <int>\n"
        "          [--outer5 <seq>] [--outer3 <seq>]\n"
        "          [--min-umi <int (default 1)>] [--max-umi <int (default 100)>]\n"
        "          [--max-errors <int (default 2)>]\n"
        "          [--threads <int (default 1)>]\n",
        prog);
    exit(1);
}

/* ── main ───────────────────────────────────────────────────────────────── */

int main(int argc, char **argv)
{
    char *input_path  = NULL;
    char *output_path = NULL;
    char *fail_path   = NULL;
    int   n_threads   = 1;
    ParseArgs pa = {
        .inner5 = NULL, .inner5_len = 0,
        .inner3 = NULL, .inner3_len = 0,
        .outer5 = NULL, .outer5_len = 0,
        .outer3 = NULL, .outer3_len = 0,
        .barcode_length = 0,
        .min_umi    = 1,
        .max_umi    = 100,
        .max_errors = 2,
    };

    static const struct option long_opts[] = {
        { "in",             required_argument, 0, 'I' },
        { "out",            required_argument, 0, 'O' },
        { "fail",           required_argument, 0, 'F' },
        { "inner5",         required_argument, 0, '5' },
        { "inner3",         required_argument, 0, '3' },
        { "outer5",         required_argument, 0, 'A' },
        { "outer3",         required_argument, 0, 'B' },
        { "barcode-length", required_argument, 0, 'L' },
        { "min-umi",        required_argument, 0, 'm' },
        { "max-umi",        required_argument, 0, 'M' },
        { "max-errors",     required_argument, 0, 'e' },
        { "threads",        required_argument, 0, 't' },
        { 0, 0, 0, 0 }
    };

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "", long_opts, &idx)) != -1) {
        switch (opt) {
            case 'I': input_path        = optarg;       break;
            case 'O': output_path       = optarg;       break;
            case 'F': fail_path         = optarg;       break;
            case '5': pa.inner5         = optarg;       break;
            case '3': pa.inner3         = optarg;       break;
            case 'A': pa.outer5         = optarg;       break;
            case 'B': pa.outer3         = optarg;       break;
            case 'L': pa.barcode_length = atoi(optarg); break;
            case 'm': pa.min_umi        = atoi(optarg); break;
            case 'M': pa.max_umi        = atoi(optarg); break;
            case 'e': pa.max_errors     = atoi(optarg); break;
            case 't': n_threads         = atoi(optarg); break;
            default:  usage(argv[0]);
        }
    }

    if (!input_path || !output_path || !fail_path ||
        !pa.inner5 || !pa.inner3 || pa.barcode_length <= 0)
        usage(argv[0]);

    if (n_threads < 1) n_threads = 1;

    pa.inner5_len = (int)strlen(pa.inner5);
    pa.inner3_len = (int)strlen(pa.inner3);
    if (pa.outer5) pa.outer5_len = (int)strlen(pa.outer5);
    if (pa.outer3) pa.outer3_len = (int)strlen(pa.outer3);

    /* ── open files ── */
    gzFile in_fp = gzopen(input_path, "r");
    if (!in_fp) { perror(input_path); return 1; }
    gzbuffer(in_fp, 1 << 20);  /* 1 MB I/O buffer */

    FILE *fail_fp = fopen(fail_path, "w");
    if (!fail_fp) { perror(fail_path); return 1; }

    /* ── initialise work queue ── */
    int queue_cap = n_threads + 2;
    if (queue_cap < MIN_QUEUE_BATCHES) queue_cap = MIN_QUEUE_BATCHES;
    if (queue_cap > MAX_QUEUE_BATCHES) queue_cap = MAX_QUEUE_BATCHES;

    WorkQueue q;
    queue_init(&q, queue_cap);

    /* ── initialise and launch worker threads ── */
    pthread_mutex_t fail_mutex;
    pthread_mutex_init(&fail_mutex, NULL);

    WorkerArgs *workers = calloc((size_t)n_threads, sizeof(WorkerArgs));
    pthread_t  *tids    = malloc((size_t)n_threads * sizeof(pthread_t));
    if (!workers || !tids) { perror("malloc"); return 1; }

    for (int i = 0; i < n_threads; i++) {
        workers[i].queue        = &q;
        workers[i].pa           = &pa;
        workers[i].local_tuples = NULL;
        workers[i].n_parsed     = 0;
        workers[i].n_failed     = 0;
        workers[i].fail_fp      = fail_fp;
        workers[i].fail_mutex   = &fail_mutex;
        if (pthread_create(&tids[i], NULL, worker_thread, &workers[i]) != 0) {
            perror("pthread_create"); return 1;
        }
    }

    /* ── producer: accumulate records into batches and push to queue ── */
    RecordBatch *cur = malloc(sizeof(RecordBatch));
    if (!cur) { perror("malloc"); return 1; }
    cur->count = 0;

    long n_total = 0;
    while (read_fastq_record(in_fp, &cur->recs[cur->count])) {
        n_total++;
        if (++cur->count == BATCH_SIZE) {
            queue_push_batch(&q, cur);
            cur->count = 0;
        }
    }
    gzclose(in_fp);

    if (cur->count > 0)          /* push any partial final batch */
        queue_push_batch(&q, cur);
    free(cur);

    queue_finish(&q);            /* signal EOF to all workers */

    /* ── join workers ── */
    for (int i = 0; i < n_threads; i++)
        pthread_join(tids[i], NULL);

    fclose(fail_fp);
    pthread_mutex_destroy(&fail_mutex);
    queue_destroy(&q);
    free(tids);

    /* ── merge per-worker hash tables ── */
    TupleEntry *tuples  = NULL;
    long n_parsed = 0, n_failed = 0;
    for (int i = 0; i < n_threads; i++) {
        merge_tables(&tuples, workers[i].local_tuples);
        n_parsed += workers[i].n_parsed;
        n_failed += workers[i].n_failed;
    }
    free(workers);

    /* ══ everything below is identical to the single-threaded version ════════ */

    int n_tuples = (int)HASH_COUNT(tuples);
    long *counts_arr = NULL;
    if (n_tuples > 0) {
        counts_arr = malloc((size_t)n_tuples * sizeof(long));
        if (!counts_arr) { perror("malloc"); return 1; }
        int i = 0;
        TupleEntry *e, *tmp;
        HASH_ITER(hh, tuples, e, tmp) counts_arr[i++] = e->count;
    }

    double g      = gini(counts_arr, n_tuples);
    long   max_ct = 0;
    double mean_ct = 0.0;
    for (int i = 0; i < n_tuples; i++) {
        if (counts_arr[i] > max_ct) max_ct = counts_arr[i];
        mean_ct += (double)counts_arr[i];
    }
    if (n_tuples > 0) mean_ct /= n_tuples;
    free(counts_arr);

    fprintf(stderr, "Threads used:           %d\n",    n_threads);
    fprintf(stderr, "Reads processed:        %ld\n",   n_total);
    fprintf(stderr, "Reads parsed:           %ld (%.1f%%)\n",
            n_parsed, n_total ? 100.0 * n_parsed / n_total : 0.0);
    fprintf(stderr, "Reads failed:           %ld (%.1f%%)\n",
            n_failed, n_total ? 100.0 * n_failed / n_total : 0.0);
    fprintf(stderr, "Unique (UMI,BC) tuples: %d\n",    n_tuples);
    fprintf(stderr, "Gini coefficient:       %.4f  (0 = no jackpotting, 1 = fully jackpotted)\n", g);
    fprintf(stderr, "Mean reads per tuple:   %.1f\n",  mean_ct);
    fprintf(stderr, "Max reads per tuple:    %ld\n",   max_ct);

    BarcodeSummary *bc_map = NULL;
    {
        TupleEntry *e, *tmp;
        HASH_ITER(hh, tuples, e, tmp) {
            const char *tab1 = strchr(e->key, '\t');
            const char *tab2 = tab1 ? strchr(tab1 + 1, '\t') : NULL;
            if (!tab1 || !tab2) continue;

            int bc_len = (int)(tab2 - tab1 - 1);
            char bc_str[MAX_BC];
            if (bc_len <= 0 || bc_len >= MAX_BC) continue;
            memcpy(bc_str, tab1 + 1, (size_t)bc_len);
            bc_str[bc_len] = '\0';

            int u5_len = (int)(tab1 - e->key);
            char umi_key[MAX_UMI_KEY];
            int uklen = snprintf(umi_key, sizeof(umi_key), "%.*s\t%s",
                                 u5_len, e->key, tab2 + 1);
            if (uklen <= 0 || uklen >= (int)sizeof(umi_key)) continue;

            BarcodeSummary *bs = NULL;
            HASH_FIND_STR(bc_map, bc_str, bs);
            if (!bs) {
                bs = calloc(1, sizeof(BarcodeSummary));
                if (!bs) { perror("calloc"); return 1; }
                memcpy(bs->barcode, bc_str, (size_t)bc_len + 1);
                HASH_ADD_STR(bc_map, barcode, bs);
            }
            bs->total_reads += e->count;

            UmiPairEntry *up = NULL;
            HASH_FIND_STR(bs->umi_set, umi_key, up);
            if (!up) {
                up = malloc(sizeof(UmiPairEntry));
                if (!up) { perror("malloc"); return 1; }
                memcpy(up->key, umi_key, (size_t)uklen + 1);
                HASH_ADD_STR(bs->umi_set, key, up);
                bs->num_unique_umis++;
            }
        }
    }

    {
        TupleEntry *e, *tmp;
        HASH_ITER(hh, tuples, e, tmp) { HASH_DEL(tuples, e); free(e); }
    }
    {
        BarcodeSummary *bs, *tmp;
        HASH_ITER(hh, bc_map, bs, tmp) {
            UmiPairEntry *up, *utmp;
            HASH_ITER(hh, bs->umi_set, up, utmp) { HASH_DEL(bs->umi_set, up); free(up); }
        }
    }

    int n_bc = (int)HASH_COUNT(bc_map);
    BarcodeSummary **bc_arr = malloc((size_t)n_bc * sizeof(BarcodeSummary *));
    if (!bc_arr && n_bc > 0) { perror("malloc"); return 1; }
    {
        int i = 0;
        BarcodeSummary *bs, *tmp;
        HASH_ITER(hh, bc_map, bs, tmp) bc_arr[i++] = bs;
    }
    qsort(bc_arr, (size_t)n_bc, sizeof(BarcodeSummary *), cmp_bc_desc);

    FILE *out_fp = fopen(output_path, "w");
    if (!out_fp) { perror(output_path); return 1; }
    fprintf(out_fp, "barcode\ttotal_reads\tnum_unique_umis\n");
    for (int i = 0; i < n_bc; i++)
        fprintf(out_fp, "%s\t%ld\t%ld\n",
                bc_arr[i]->barcode,
                bc_arr[i]->total_reads,
                bc_arr[i]->num_unique_umis);
    fclose(out_fp);
    free(bc_arr);

    fprintf(stderr, "Unique barcodes:         %d\n",  n_bc);
    fprintf(stderr, "Results written to:      %s\n",  output_path);
    fprintf(stderr, "Failed reads written to: %s\n",  fail_path);

    {
        BarcodeSummary *bs, *tmp;
        HASH_ITER(hh, bc_map, bs, tmp) { HASH_DEL(bc_map, bs); free(bs); }
    }

    return 0;
}
