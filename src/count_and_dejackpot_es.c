/*
 * count_and_dejackpot_es.c  –  memory-bounded external-sort variant
 *
 * Identical logic to count_and_dejackpot_mt.c for the parse phase, but replaces
 * the in-memory hash table with an external sort + k-way merge approach so that
 * peak RAM usage is proportional to --sort-mem regardless of dataset size.
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
 *   2. pthreads   – POSIX threads (included with glibc; link with -lpthread)
 *
 *   Note: uthash.h is NOT required for this variant.
 *
 * Recommended build commands:
 *   gcc  -std=c99 -O2 -o count_and_dejackpot_es count_and_dejackpot_es.c -lz -lm -lpthread
 *   clang -std=c99 -O2 -o count_and_dejackpot_es count_and_dejackpot_es.c -lz -lm -lpthread
 *
 * Note: -std=c99 is required.  Older GCC installations (common on clusters)
 * default to C89.  Use -std=gnu99 if you need GNU extensions (strdup, getline).
 *
 * ── ALGORITHM ─────────────────────────────────────────────────────────────────
 *
 * Phase 1 – Parse + sort runs (multi-threaded):
 *   Workers parse FASTQ records and accumulate sort keys into a per-worker
 *   flat buffer.  Each key is "bc\tumi5\tumi3" (barcode-primary).  When a
 *   worker's buffer fills (sort_mem / n_threads bytes), it sorts in place and
 *   flushes to a numbered temp file.  Workers continue into a fresh buffer.
 *
 * Phase 2 – k-way merge (serial, main thread):
 *   All temp files are opened and their heads loaded into a min-heap.  The
 *   heap produces a globally sorted stream of keys.  Consecutive identical
 *   keys are the same (bc, umi5, umi3) tuple seen in multiple reads; runs of
 *   identical keys are collapsed into a per-tuple count.  Per-barcode totals
 *   and unique-UMI counts are accumulated directly from the sorted stream;
 *   no second hash table is needed.
 *
 * Phase 3 – Output:
 *   Barcode rows are collected during the merge, then sorted by num_unique_umis
 *   descending before writing to the output TSV.
 *
 * ── MEMORY ────────────────────────────────────────────────────────────────────
 *
 *   Sort buffers:  --sort-mem MB total (split evenly across workers)
 *   Gini sample :  --gini-cap × 8 bytes (default 8 MB)
 *   Barcode rows:  O(num_unique_barcodes × 272 bytes)  — usually small
 *   Merge heads  :  O(n_temp_files × ~5 KB line buffer)  — very small
 *
 * ── THREADING DESIGN ──────────────────────────────────────────────────────────
 *
 *   Same producer / consumer queue as count_and_dejackpot_mt.c.
 *   Workers now fill a private sort buffer instead of a hash table.
 *   The merge phase is serial.
 *
 * ── FUZZY MATCHING NOTE ───────────────────────────────────────────────────────
 *
 *   Substitutions only (Hamming distance).  See count_and_dejackpot.c for
 *   rationale.
 */

/* Required for mkdtemp, strdup, PATH_MAX, sysconf on POSIX systems */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <pthread.h>
#include <zlib.h>
#include <unistd.h>   /* sysconf */
#include <sys/stat.h> /* stat   */
#include <limits.h>   /* PATH_MAX */

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

/* ── constants ──────────────────────────────────────────────────────────── */

#define MAX_LINE    16384  /* max bytes in a single FASTQ line              */
#define MAX_UMI       256  /* max UMI sequence length                       */
#define MAX_BC        256  /* max barcode sequence length                   */

/* Sort key format: "bc\tumi5\tumi3\0"  (barcode-primary for sorted merge) */
#define MAX_SORT_KEY  (MAX_BC + 1 + MAX_UMI + 1 + MAX_UMI + 1)

#define BATCH_SIZE          32     /* FASTQ records per queue slot           */
#define MIN_QUEUE_BATCHES    4     /* minimum queue depth (in batches)       */
#define MAX_QUEUE_BATCHES   64     /* maximum queue depth (in batches)       */
#define FAIL_BUF_BYTES   (256 * 1024)  /* per-thread fail-record buffer     */

#define SORT_MEM_DEFAULT   4096L  /* total sort memory, MB                  */
#define GINI_CAP_DEFAULT   1000000L  /* max tuples sampled for Gini         */

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

typedef struct {
    char umi5   [MAX_UMI]; int u5len;
    char barcode[MAX_BC];  int bclen;
    char umi3   [MAX_UMI]; int u3len;
} ParseResult;

/* ── FASTQ record and batch ─────────────────────────────────────────────── */

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
    int             capacity;
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

static inline void copy_record(FastqRecord *dst, const FastqRecord *src)
{
    dst->hlen = src->hlen; memcpy(dst->header, src->header, (size_t)src->hlen + 1);
    dst->slen = src->slen; memcpy(dst->seq,    src->seq,    (size_t)src->slen + 1);
    dst->qlen = src->qlen; memcpy(dst->qual,   src->qual,   (size_t)src->qlen + 1);
}

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

static void fail_buf_write(FailBuffer *fb,
                            const char *hdr, int hlen,
                            const char *seq, int slen,
                            const char *qual, int qlen,
                            FILE *fp, pthread_mutex_t *mu)
{
    int need = hlen + slen + qlen + 5;
    if (fb->used + need > FAIL_BUF_BYTES)
        fail_buf_flush(fb, fp, mu);
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

/* ── sort state (per-worker sort buffer + temp file list) ───────────────── */

/*
 * Each entry in the sort buffer is exactly `stride` bytes long.
 * The first bytes hold a null-terminated sort key "bc\tumi5\tumi3".
 * Bytes beyond the null terminator are unused padding; qsort passes
 * element-start pointers to cmp_sort_entry which just calls strcmp.
 */
typedef struct {
    char    *buf;        /* flat array: n × stride bytes                   */
    long     n;          /* entries currently in buf                       */
    long     cap;        /* max entries (= mem_bytes / stride)             */
    int      stride;     /* bytes per entry                                */
    char   **paths;      /* paths of temp files written by this worker     */
    int      n_paths;
    int      paths_cap;
    char     tmp_dir[PATH_MAX];
    int      worker_id;
} SortState;

/* Comparator for qsort: each element is stride bytes; key starts at byte 0 */
static int cmp_sort_entry(const void *a, const void *b)
{
    return strcmp((const char *)a, (const char *)b);
}

static void sort_state_flush(SortState *ss)
{
    if (ss->n == 0) return;

    qsort(ss->buf, (size_t)ss->n, (size_t)ss->stride, cmp_sort_entry);

    /* PATH_MAX is enough for tmp_dir + "/run_NNNN_NNNN.tmp" (19 chars) */
    char path[PATH_MAX + 32];
    snprintf(path, sizeof(path), "%s/run_%04d_%04d.tmp",
             ss->tmp_dir, ss->worker_id, ss->n_paths);
    FILE *fp = fopen(path, "w");
    if (!fp) { perror(path); exit(1); }

    for (long i = 0; i < ss->n; i++) {
        fputs(ss->buf + i * ss->stride, fp);
        fputc('\n', fp);
    }
    fclose(fp);

    /* Record path for later merge */
    if (ss->n_paths == ss->paths_cap) {
        ss->paths_cap *= 2;
        ss->paths = realloc(ss->paths, (size_t)ss->paths_cap * sizeof(char *));
        if (!ss->paths) { perror("realloc paths"); exit(1); }
    }
    ss->paths[ss->n_paths] = strdup(path);
    if (!ss->paths[ss->n_paths]) { perror("strdup"); exit(1); }
    ss->n_paths++;
    ss->n = 0;
}

static void sort_state_init(SortState *ss, long mem_bytes, int stride,
                             const char *tmp_dir, int worker_id)
{
    ss->stride = stride;
    ss->cap    = mem_bytes / stride;
    if (ss->cap < 1) ss->cap = 1;

    ss->buf = malloc((size_t)ss->cap * (size_t)stride);
    if (!ss->buf) { perror("malloc sort buffer"); exit(1); }

    ss->n          = 0;
    ss->paths_cap  = 64;
    ss->n_paths    = 0;
    ss->paths = malloc((size_t)ss->paths_cap * sizeof(char *));
    if (!ss->paths) { perror("malloc paths"); exit(1); }

    snprintf(ss->tmp_dir, sizeof(ss->tmp_dir), "%s", tmp_dir);
    ss->worker_id = worker_id;
}

/* Called from worker after popping last batch to flush any remaining tuples */
static void sort_state_finish(SortState *ss)
{
    sort_state_flush(ss);
}

static void sort_state_free(SortState *ss)
{
    free(ss->buf);
    for (int i = 0; i < ss->n_paths; i++) free(ss->paths[i]);
    free(ss->paths);
}

/* Add one parsed result to the sort buffer; flushes to a temp file if full */
static void sort_state_add(SortState *ss, const ParseResult *pr)
{
    if (ss->n == ss->cap) sort_state_flush(ss);

    /* Build sort key: bc\tumi5\tumi3\0 */
    char *entry = ss->buf + ss->n * ss->stride;
    char *p = entry;
    memcpy(p, pr->barcode, (size_t)pr->bclen); p += pr->bclen; *p++ = '\t';
    memcpy(p, pr->umi5,    (size_t)pr->u5len); p += pr->u5len; *p++ = '\t';
    memcpy(p, pr->umi3,    (size_t)pr->u3len); p += pr->u3len; *p   = '\0';
    ss->n++;
}

/* ── worker thread state ────────────────────────────────────────────────── */

typedef struct {
    WorkQueue       *queue;
    const ParseArgs *pa;
    SortState       *ss;
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

    memcpy(out->umi5,    seq + u5_start,          (size_t)u5_len);
    out->umi5[u5_len]    = '\0';
    memcpy(out->barcode, seq + i5_end,  (size_t)pa->barcode_length);
    out->barcode[pa->barcode_length] = '\0';
    memcpy(out->umi3,    seq + u3_start,          (size_t)u3_len);
    out->umi3[u3_len]    = '\0';
    out->u5len = u5_len;
    out->bclen = pa->barcode_length;
    out->u3len = u3_len;
    return 0;
}

/* ── worker thread ──────────────────────────────────────────────────────── */

static void *worker_thread(void *arg)
{
    WorkerArgs *wa = arg;

    RecordBatch *batch = malloc(sizeof(RecordBatch));
    FailBuffer   fb    = { malloc(FAIL_BUF_BYTES), 0 };
    if (!batch || !fb.data) { perror("malloc"); exit(1); }

    while (queue_pop_batch(wa->queue, batch)) {
        for (int i = 0; i < batch->count; i++) {
            const FastqRecord *r = &batch->recs[i];
            ParseResult pr;

            if (parse_read(r->seq, wa->pa, &pr) == 0) {
                sort_state_add(wa->ss, &pr);
                wa->n_parsed++;
            } else {
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
    sort_state_finish(wa->ss);  /* flush any remaining tuples */
    free(fb.data);
    free(batch);
    return NULL;
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

/* ── merge stream + min-heap ────────────────────────────────────────────── */

typedef struct {
    FILE *fp;
    char  cur[MAX_SORT_KEY + 2];  /* current line, newline stripped */
    int   cur_len;
    char  path[PATH_MAX];
    int   exhausted;
} MergeStream;

/* Read next line from stream; returns 1 on success, 0 on EOF */
static int stream_advance(MergeStream *ms)
{
    if (!fgets(ms->cur, (int)sizeof(ms->cur), ms->fp)) {
        ms->exhausted = 1;
        ms->cur[0]    = '\0';
        ms->cur_len   = 0;
        return 0;
    }
    ms->cur_len = (int)strlen(ms->cur);
    while (ms->cur_len > 0 &&
           (ms->cur[ms->cur_len - 1] == '\n' || ms->cur[ms->cur_len - 1] == '\r'))
        ms->cur[--ms->cur_len] = '\0';
    return 1;
}

typedef struct {
    MergeStream **arr;
    int           n;
    int           cap;
} MinHeap;

static void heap_push(MinHeap *h, MergeStream *ms)
{
    if (h->n == h->cap) {
        h->cap = h->cap ? h->cap * 2 : 16;
        h->arr = realloc(h->arr, (size_t)h->cap * sizeof(MergeStream *));
        if (!h->arr) { perror("realloc heap"); exit(1); }
    }
    int i = h->n++;
    h->arr[i] = ms;
    /* Sift up */
    while (i > 0) {
        int parent = (i - 1) / 2;
        if (strcmp(h->arr[parent]->cur, h->arr[i]->cur) <= 0) break;
        MergeStream *tmp  = h->arr[parent];
        h->arr[parent]    = h->arr[i];
        h->arr[i]         = tmp;
        i = parent;
    }
}

static MergeStream *heap_pop(MinHeap *h)
{
    MergeStream *ret = h->arr[0];
    h->arr[0] = h->arr[--h->n];
    /* Sift down */
    int i = 0;
    for (;;) {
        int left = 2*i+1, right = 2*i+2, smallest = i;
        if (left  < h->n && strcmp(h->arr[left]->cur,    h->arr[smallest]->cur) < 0) smallest = left;
        if (right < h->n && strcmp(h->arr[right]->cur, h->arr[smallest]->cur) < 0) smallest = right;
        if (smallest == i) break;
        MergeStream *tmp    = h->arr[i];
        h->arr[i]           = h->arr[smallest];
        h->arr[smallest]    = tmp;
        i = smallest;
    }
    return ret;
}

/* ── per-barcode summary (collected during merge) ───────────────────────── */

typedef struct {
    char barcode[MAX_BC];
    long total_reads;
    long num_unique_umis;
} BcSummary;

static int cmp_bc_desc(const void *a, const void *b)
{
    const BcSummary *x = (const BcSummary *)a;
    const BcSummary *y = (const BcSummary *)b;
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
        "          [--threads <int (default 1)>]\n"
        "          [--sort-mem <MB (default 4096)>]\n"
        "          [--gini-cap <N tuples (default 1000000)>]\n"
        "          [--tmp-dir <path (default '.')>]\n",
        prog);
    exit(1);
}

/* ── main ───────────────────────────────────────────────────────────────── */

int main(int argc, char **argv)
{
    char *input_path    = NULL;
    char *output_path   = NULL;
    char *fail_path     = NULL;
    char *tmp_dir_base  = ".";
    int   n_threads     = 1;
    long  sort_mem_mb   = SORT_MEM_DEFAULT;
    long  gini_cap      = GINI_CAP_DEFAULT;

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
        { "sort-mem",       required_argument, 0, 's' },
        { "gini-cap",       required_argument, 0, 'g' },
        { "tmp-dir",        required_argument, 0, 'd' },
        { 0, 0, 0, 0 }
    };

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "", long_opts, &idx)) != -1) {
        switch (opt) {
            case 'I': input_path        = optarg;           break;
            case 'O': output_path       = optarg;           break;
            case 'F': fail_path         = optarg;           break;
            case '5': pa.inner5         = optarg;           break;
            case '3': pa.inner3         = optarg;           break;
            case 'A': pa.outer5         = optarg;           break;
            case 'B': pa.outer3         = optarg;           break;
            case 'L': pa.barcode_length = atoi(optarg);     break;
            case 'm': pa.min_umi        = atoi(optarg);     break;
            case 'M': pa.max_umi        = atoi(optarg);     break;
            case 'e': pa.max_errors     = atoi(optarg);     break;
            case 't': n_threads         = atoi(optarg);     break;
            case 's': sort_mem_mb       = atol(optarg);     break;
            case 'g': gini_cap          = atol(optarg);     break;
            case 'd': tmp_dir_base      = optarg;           break;
            default:  usage(argv[0]);
        }
    }

    if (!input_path || !output_path || !fail_path ||
        !pa.inner5 || !pa.inner3 || pa.barcode_length <= 0)
        usage(argv[0]);

    if (n_threads  < 1) n_threads  = 1;
    if (sort_mem_mb < 1) sort_mem_mb = 1;
    if (gini_cap   < 1) gini_cap   = 1;

    pa.inner5_len = (int)strlen(pa.inner5);
    pa.inner3_len = (int)strlen(pa.inner3);
    if (pa.outer5) pa.outer5_len = (int)strlen(pa.outer5);
    if (pa.outer3) pa.outer3_len = (int)strlen(pa.outer3);

    /*
     * stride = worst-case sort key length + 1.
     * Sort key: "bc\tumi5\tumi3\0"
     * = barcode_length + 1 + max_umi + 1 + max_umi + 1
     *                        (tab)          (tab)    (null)
     */
    int stride = pa.barcode_length + 1 + pa.max_umi + 1 + pa.max_umi + 1;
    if (stride >= MAX_SORT_KEY) {
        fprintf(stderr, "Error: barcode_length + 2*max_umi too large "
                "(stride=%d, max=%d).\n", stride, MAX_SORT_KEY - 1);
        return 1;
    }

    long sort_mem_bytes      = sort_mem_mb * 1024L * 1024L;
    long sort_mem_per_worker = sort_mem_bytes / n_threads;
    if (sort_mem_per_worker < stride) {
        fprintf(stderr,
            "Warning: sort buffer per worker (%ld bytes) is smaller than one "
            "sort key (%d bytes).\nConsider increasing --sort-mem.\n",
            sort_mem_per_worker, stride);
        sort_mem_per_worker = stride;
    }

    /* ── create temp directory ── */
    char tmp_dir[PATH_MAX];
    snprintf(tmp_dir, sizeof(tmp_dir), "%s/dejackpot_tmp_XXXXXX", tmp_dir_base);
    if (!mkdtemp(tmp_dir)) {
        fprintf(stderr, "Error: could not create temp directory in '%s': ", tmp_dir_base);
        perror("");
        return 1;
    }
    fprintf(stderr, "Temp directory: %s\n", tmp_dir);

    /* ── open input / fail files ── */
    gzFile in_fp = gzopen(input_path, "r");
    if (!in_fp) { perror(input_path); rmdir(tmp_dir); return 1; }
    gzbuffer(in_fp, 1 << 20);  /* 1 MB I/O buffer */

    FILE *fail_fp = fopen(fail_path, "w");
    if (!fail_fp) { perror(fail_path); gzclose(in_fp); rmdir(tmp_dir); return 1; }

    /* ── initialise work queue ── */
    int queue_cap = n_threads + 2;
    if (queue_cap < MIN_QUEUE_BATCHES) queue_cap = MIN_QUEUE_BATCHES;
    if (queue_cap > MAX_QUEUE_BATCHES) queue_cap = MAX_QUEUE_BATCHES;

    WorkQueue q;
    queue_init(&q, queue_cap);

    /* ── initialise per-worker sort state ── */
    SortState *sort_states = malloc((size_t)n_threads * sizeof(SortState));
    if (!sort_states) { perror("malloc sort_states"); return 1; }
    for (int i = 0; i < n_threads; i++)
        sort_state_init(&sort_states[i], sort_mem_per_worker, stride, tmp_dir, i);

    /* ── initialise and launch worker threads ── */
    pthread_mutex_t fail_mutex;
    pthread_mutex_init(&fail_mutex, NULL);

    WorkerArgs *workers = calloc((size_t)n_threads, sizeof(WorkerArgs));
    pthread_t  *tids    = malloc((size_t)n_threads * sizeof(pthread_t));
    if (!workers || !tids) { perror("malloc"); return 1; }

    for (int i = 0; i < n_threads; i++) {
        workers[i].queue       = &q;
        workers[i].pa          = &pa;
        workers[i].ss          = &sort_states[i];
        workers[i].n_parsed    = 0;
        workers[i].n_failed    = 0;
        workers[i].fail_fp     = fail_fp;
        workers[i].fail_mutex  = &fail_mutex;
        if (pthread_create(&tids[i], NULL, worker_thread, &workers[i]) != 0) {
            perror("pthread_create"); return 1;
        }
    }

    /* ── producer: read FASTQ and push batches ── */
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

    if (cur->count > 0)
        queue_push_batch(&q, cur);
    free(cur);

    queue_finish(&q);

    /* ── join workers (workers flush final sort runs here) ── */
    for (int i = 0; i < n_threads; i++)
        pthread_join(tids[i], NULL);

    fclose(fail_fp);
    pthread_mutex_destroy(&fail_mutex);
    queue_destroy(&q);
    free(tids);

    /* ── collect stats and temp file list ── */
    long n_parsed = 0, n_failed = 0;
    int  n_total_runs = 0;
    for (int i = 0; i < n_threads; i++) {
        n_parsed     += workers[i].n_parsed;
        n_failed     += workers[i].n_failed;
        n_total_runs += sort_states[i].n_paths;
    }
    free(workers);

    /* ── fd budget check ── */
    {
        long fd_limit    = sysconf(_SC_OPEN_MAX);
        long fd_reserved = 20;  /* stdin/stdout/stderr/gzip/output/fail/etc. */
        if (n_total_runs > fd_limit - fd_reserved) {
            fprintf(stderr,
                "Error: %d sort runs would require %d file descriptors, but only "
                "%ld are available (fd limit: %ld).\n"
                "  Increase --sort-mem (current: %ld MB) or run: ulimit -n %d\n",
                n_total_runs, n_total_runs,
                fd_limit - fd_reserved, fd_limit,
                sort_mem_mb, n_total_runs + (int)fd_reserved);
            /* Attempt cleanup */
            for (int i = 0; i < n_threads; i++) {
                for (int j = 0; j < sort_states[i].n_paths; j++)
                    remove(sort_states[i].paths[j]);
                sort_state_free(&sort_states[i]);
            }
            free(sort_states);
            rmdir(tmp_dir);
            return 1;
        }
    }

    /* ── handle empty dataset ── */
    if (n_total_runs == 0) {
        fprintf(stderr, "Threads used:           %d\n", n_threads);
        fprintf(stderr, "Reads processed:        %ld\n", n_total);
        fprintf(stderr, "Reads parsed:           0 (0.0%%)\n");
        fprintf(stderr, "Reads failed:           %ld (%.1f%%)\n",
                n_failed, n_total ? 100.0 * n_failed / n_total : 0.0);
        fprintf(stderr, "Unique (UMI,BC) tuples: 0\n");
        fprintf(stderr, "Gini coefficient:       nan\n");
        fprintf(stderr, "Unique barcodes:         0\n");

        FILE *out_fp = fopen(output_path, "w");
        if (out_fp) {
            fprintf(out_fp, "barcode\ttotal_reads\tnum_unique_umis\n");
            fclose(out_fp);
        }
        for (int i = 0; i < n_threads; i++) sort_state_free(&sort_states[i]);
        free(sort_states);
        rmdir(tmp_dir);
        return 0;
    }

    /* ── open all temp files and build initial heap ── */
    MergeStream *streams = malloc((size_t)n_total_runs * sizeof(MergeStream));
    if (!streams) { perror("malloc streams"); return 1; }

    int k = 0;
    for (int i = 0; i < n_threads; i++) {
        for (int j = 0; j < sort_states[i].n_paths; j++) {
            const char *path = sort_states[i].paths[j];
            streams[k].fp = fopen(path, "r");
            if (!streams[k].fp) {
                perror(path);
                /* Cleanup already-opened streams */
                for (int m = 0; m < k; m++) fclose(streams[m].fp);
                free(streams);
                for (int m = 0; m < n_threads; m++) {
                    for (int n = 0; n < sort_states[m].n_paths; n++)
                        remove(sort_states[m].paths[n]);
                    sort_state_free(&sort_states[m]);
                }
                free(sort_states);
                rmdir(tmp_dir);
                return 1;
            }
            strncpy(streams[k].path, path, sizeof(streams[k].path) - 1);
            streams[k].path[sizeof(streams[k].path) - 1] = '\0';
            streams[k].exhausted = 0;
            stream_advance(&streams[k]);
            k++;
        }
    }

    /* Free sort_states (temp file paths are now in streams[].path) */
    for (int i = 0; i < n_threads; i++) sort_state_free(&sort_states[i]);
    free(sort_states);

    MinHeap heap = { NULL, 0, 0 };
    for (int i = 0; i < n_total_runs; i++) {
        if (!streams[i].exhausted)
            heap_push(&heap, &streams[i]);
    }

    /* ── allocate Gini sample array and barcode summary array ── */
    long  *gini_counts = malloc((size_t)gini_cap * sizeof(long));
    int    gini_n      = 0;
    if (!gini_counts) { perror("malloc gini_counts"); return 1; }

    int       bc_cap = 1024;
    BcSummary *bc_arr = malloc((size_t)bc_cap * sizeof(BcSummary));
    int        bc_n   = 0;
    if (!bc_arr) { perror("malloc bc_arr"); return 1; }

    /* ── k-way merge loop ────────────────────────────────────────────────
     *
     * State machine:
     *   prev_key     – sort key of the last popped tuple
     *   prev_bc      – barcode of the last processed barcode
     *   prev_bc_len  – strlen(prev_bc) cached to avoid repeated strlen
     *   cur_tuple_count – consecutive reads with the current tuple key
     *   cur_total    – total reads seen for the current barcode
     *   cur_unique   – unique (umi5, umi3) pairs for the current barcode
     *
     * Transition rules:
     *   key == prev_key: same tuple, cur_tuple_count++, cur_total++
     *   key != prev_key:
     *     emit cur_tuple_count to Gini sample (if under cap)
     *     same barcode → cur_total++, cur_unique++
     *     new barcode  → emit barcode row, reset counters
     * ──────────────────────────────────────────────────────────────────── */

    char prev_key[MAX_SORT_KEY] = "";
    int  prev_key_len           = 0;
    char prev_bc [MAX_BC]       = "";
    int  prev_bc_len            = 0;
    long cur_tuple_count        = 0;
    long cur_total              = 0;
    long cur_unique             = 0;
    long n_tuples_total         = 0;
    long max_ct                 = 0;
    long total_ct_sum           = 0;   /* sum of all tuple counts = n_parsed */

    while (heap.n > 0) {
        MergeStream *ms  = heap_pop(&heap);
        const char  *key = ms->cur;
        int          klen = ms->cur_len;

        /* Extract barcode (prefix up to first '\t') */
        const char *tab1  = memchr(key, '\t', (size_t)klen);
        int         bclen = tab1 ? (int)(tab1 - key) : klen;

        if (klen == prev_key_len && memcmp(key, prev_key, (size_t)klen) == 0) {
            /* ── same tuple as previous read ── */
            cur_tuple_count++;
            cur_total++;
        } else {
            /* ── new tuple ── */
            if (cur_tuple_count > 0) {
                /* Emit previous tuple count */
                if (gini_n < (int)gini_cap)
                    gini_counts[gini_n++] = cur_tuple_count;
                n_tuples_total++;
                if (cur_tuple_count > max_ct) max_ct = cur_tuple_count;
                total_ct_sum += cur_tuple_count;
            }
            cur_tuple_count = 1;

            /* Check whether barcode changed */
            int bc_changed = (bclen != prev_bc_len ||
                              memcmp(key, prev_bc, (size_t)bclen) != 0);

            if (bc_changed) {
                if (prev_bc_len > 0) {
                    /* Emit completed barcode row */
                    if (bc_n == bc_cap) {
                        bc_cap *= 2;
                        bc_arr = realloc(bc_arr, (size_t)bc_cap * sizeof(BcSummary));
                        if (!bc_arr) { perror("realloc bc_arr"); return 1; }
                    }
                    memcpy(bc_arr[bc_n].barcode, prev_bc, (size_t)prev_bc_len + 1);
                    bc_arr[bc_n].total_reads    = cur_total;
                    bc_arr[bc_n].num_unique_umis = cur_unique;
                    bc_n++;
                }
                /* Reset for new barcode */
                memcpy(prev_bc, key, (size_t)bclen);
                prev_bc[bclen] = '\0';
                prev_bc_len    = bclen;
                cur_total      = 1;
                cur_unique     = 1;
            } else {
                /* Same barcode, new (umi5, umi3) combination */
                cur_total++;
                cur_unique++;
            }

            /* Update prev_key */
            memcpy(prev_key, key, (size_t)klen + 1);
            prev_key_len = klen;
        }

        /* Advance the stream; re-insert into heap or discard if exhausted */
        if (stream_advance(ms)) {
            heap_push(&heap, ms);
        } else {
            /* Stream exhausted — close and delete temp file */
            fclose(ms->fp);
            remove(ms->path);
        }
    }

    /* ── flush final tuple and barcode ── */
    if (cur_tuple_count > 0) {
        if (gini_n < (int)gini_cap) gini_counts[gini_n++] = cur_tuple_count;
        n_tuples_total++;
        if (cur_tuple_count > max_ct) max_ct = cur_tuple_count;
        total_ct_sum += cur_tuple_count;
    }
    if (prev_bc_len > 0) {
        if (bc_n == bc_cap) {
            bc_cap *= 2;
            bc_arr = realloc(bc_arr, (size_t)bc_cap * sizeof(BcSummary));
            if (!bc_arr) { perror("realloc bc_arr"); return 1; }
        }
        memcpy(bc_arr[bc_n].barcode, prev_bc, (size_t)prev_bc_len + 1);
        bc_arr[bc_n].total_reads     = cur_total;
        bc_arr[bc_n].num_unique_umis = cur_unique;
        bc_n++;
    }

    /* ── compute statistics ── */
    double g        = gini(gini_counts, gini_n);
    double mean_ct  = n_tuples_total > 0 ? (double)total_ct_sum / n_tuples_total : 0.0;
    free(gini_counts);

    fprintf(stderr, "Threads used:           %d\n",  n_threads);
    fprintf(stderr, "Reads processed:        %ld\n", n_total);
    fprintf(stderr, "Reads parsed:           %ld (%.1f%%)\n",
            n_parsed, n_total ? 100.0 * n_parsed / n_total : 0.0);
    fprintf(stderr, "Reads failed:           %ld (%.1f%%)\n",
            n_failed, n_total ? 100.0 * n_failed / n_total : 0.0);
    fprintf(stderr, "Unique (UMI,BC) tuples: %ld\n", n_tuples_total);
    if (gini_n < (int)n_tuples_total)
        fprintf(stderr, "Gini coefficient:       %.4f  (estimated from first %d of %ld tuples; "
                "0 = no jackpotting, 1 = fully jackpotted)\n", g, gini_n, n_tuples_total);
    else
        fprintf(stderr, "Gini coefficient:       %.4f  (0 = no jackpotting, 1 = fully jackpotted)\n", g);
    fprintf(stderr, "Mean reads per tuple:   %.1f\n", mean_ct);
    fprintf(stderr, "Max reads per tuple:    %ld\n",  max_ct);

    /* ── sort barcodes by unique UMIs (descending) and write output ── */
    qsort(bc_arr, (size_t)bc_n, sizeof(BcSummary), cmp_bc_desc);

    FILE *out_fp = fopen(output_path, "w");
    if (!out_fp) { perror(output_path); return 1; }
    fprintf(out_fp, "barcode\ttotal_reads\tnum_unique_umis\n");
    for (int i = 0; i < bc_n; i++)
        fprintf(out_fp, "%s\t%ld\t%ld\n",
                bc_arr[i].barcode,
                bc_arr[i].total_reads,
                bc_arr[i].num_unique_umis);
    fclose(out_fp);

    fprintf(stderr, "Unique barcodes:         %d\n", bc_n);
    fprintf(stderr, "Results written to:      %s\n", output_path);
    fprintf(stderr, "Failed reads written to: %s\n", fail_path);

    /* ── cleanup ── */
    free(bc_arr);
    free(heap.arr);
    free(streams);
    rmdir(tmp_dir);  /* should be empty; temp files removed as streams exhausted */

    return 0;
}
