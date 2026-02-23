/*
 * count_and_dejackpot.c  –  C port of count_and_dejackpot.py
 *
 * Counts barcodes in FASTQ data using pseudoUMIs to detect PCR jackpotting.
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
 * Recommended build commands:
 *   gcc  -std=c99 -O2 -o count_and_dejackpot count_and_dejackpot.c -lz -lm
 *   clang -std=c99 -O2 -o count_and_dejackpot count_and_dejackpot.c -lz -lm
 *
 * Note: -std=c99 is required on older GCC installations (common on clusters)
 * that default to C89.  Use -std=gnu99 if you need GNU extensions.
 *
 * If uthash.h is in a non-standard location add -I/path/to/uthash/.
 * For debug builds add -g -fsanitize=address,undefined.
 *
 * ── THREADING NOTE ───────────────────────────────────────────────────────────
 *
 * This version is intentionally single-threaded.  The per-read logic is isolated
 * in parse_read(), which is:
 *   • stateless  – touches no global variables
 *   • thread-safe – all working memory lives on the stack or in *out
 *
 * The read loop in main() is structured so that parallelisation requires only:
 *   1. A shared work queue fed by a reader thread (header/seq/qual records).
 *   2. N worker threads that each call parse_read() and accumulate into a
 *      *thread-local* TupleEntry hash table.
 *   3. A serial merge step that sums all local tables into one global table.
 * No other restructuring is needed; every function below is already safe to
 * call from multiple threads simultaneously on distinct data.
 *
 * ── FUZZY MATCHING NOTE ──────────────────────────────────────────────────────
 *
 * The Python version uses the `regex` library with {e<=N} which allows
 * substitutions, insertions, and deletions.  This C version allows substitutions
 * only (Hamming distance).  For DNA sequencing data, errors in constant sequences
 * are overwhelmingly substitutions, so this is sufficient in practice.
 *
 * find_core() jointly minimises the total substitution count across inner5 and
 * inner3, which matches the BESTMATCH semantics of the Python regex.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <zlib.h>
#include "uthash.h"

/* ── constants ──────────────────────────────────────────────────────────── */

#define MAX_LINE    16384  /* max bytes in a single FASTQ line (incl. newline) */
#define MAX_PAT       512  /* max constant-sequence length                     */
#define MAX_UMI       256  /* max UMI sequence length                          */
#define MAX_BC        256  /* max barcode sequence length                      */

/* key layouts stored inside hash structs */
#define MAX_TUPLE_KEY  (MAX_UMI + 1 + MAX_BC + 1 + MAX_UMI + 1)  /* umi5\tbc\tumi3\0  */
#define MAX_UMI_KEY    (MAX_UMI + 1 + MAX_UMI + 1)               /* umi5\tumi3\0      */

/* ── data structures ────────────────────────────────────────────────────── */

/* All parameters needed for per-read parsing.
 * Passed as a single const pointer so worker threads share one read-only copy. */
typedef struct {
    const char *inner5;  int inner5_len;
    const char *inner3;  int inner3_len;
    const char *outer5;  int outer5_len;  /* NULL when not supplied */
    const char *outer3;  int outer3_len;  /* NULL when not supplied */
    int barcode_length;
    int min_umi, max_umi;
    int max_errors;
} ParseArgs;

/* Output of a successful parse_read() call. */
typedef struct {
    char umi5   [MAX_UMI];
    char barcode[MAX_BC];
    char umi3   [MAX_UMI];
} ParseResult;

/* (umi5, barcode, umi3) → read count */
typedef struct TupleEntry {
    char           key[MAX_TUPLE_KEY];
    long           count;
    UT_hash_handle hh;
} TupleEntry;

/* Set element: a unique (umi5, umi3) pair within one barcode */
typedef struct UmiPairEntry {
    char           key[MAX_UMI_KEY];
    UT_hash_handle hh;
} UmiPairEntry;

/* Per-barcode summary (built in the second pass) */
typedef struct BarcodeSummary {
    char           barcode[MAX_BC];
    long           total_reads;
    long           num_unique_umis;
    UmiPairEntry  *umi_set;   /* hash-set used to deduplicate UMI pairs; freed after counting */
    UT_hash_handle hh;
} BarcodeSummary;

/* ── inner helpers: substitution counting ───────────────────────────────── */

/*
 * count_subst: count mismatches between text[pos..pos+plen) and pattern.
 * Stops early and returns limit+1 once the budget is exceeded.
 * Inline so the hot inner loop in find_core() can be optimised aggressively.
 */
static inline int count_subst(const char *text, int pos,
                               const char *pattern, int plen, int limit)
{
    int err = 0;
    for (int j = 0; j < plen; j++)
        if (text[pos + j] != pattern[j])
            if (++err > limit) return err;
    return err;
}

/* ── core-pattern search (joint BESTMATCH over inner5 + inner3) ─────────── */

/*
 * find_core: slide the composite pattern  [inner5][barcode_length chars][inner3]
 * along seq[0..slen) and find the position where the total substitutions in
 * inner5 and inner3 combined is minimised, subject to each being <= max_errors.
 *
 * On success fills *i5_start, *i5_end, *i3_start, *i3_end and returns 0.
 * Returns -1 if no valid position exists.
 *
 * Thread-safe: stateless, stack-only.
 */
static int find_core(const char *seq, int slen, const ParseArgs *pa,
                     int *i5_start, int *i5_end,
                     int *i3_start, int *i3_end)
{
    int span = pa->inner5_len + pa->barcode_length + pa->inner3_len;
    if (span > slen) return -1;

    int best_pos = -1;
    int best_err = 2 * pa->max_errors + 1;  /* sentinel: worse than any valid match */

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
            if (total == 0) break;  /* perfect match – can't improve */
        }
    }

    if (best_pos < 0) return -1;

    *i5_start = best_pos;
    *i5_end   = best_pos + pa->inner5_len;
    *i3_start = *i5_end  + pa->barcode_length;
    *i3_end   = *i3_start + pa->inner3_len;
    return 0;
}

/* ── fuzzy search for a single pattern ─────────────────────────────────── */

/*
 * fuzzy_find: find pattern[0..plen) in text[0..tlen) allowing up to
 * max_errors substitutions.  Returns the start position of the best
 * (leftmost, fewest errors) match, or -1 on failure.
 * Sets *match_end = start + plen when non-NULL and a match is found.
 *
 * Thread-safe: stateless, stack-only.
 */
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

/* ── read parsing ───────────────────────────────────────────────────────── */

/*
 * parse_read: extract umi5, barcode, umi3 from `seq` using the parameters in `pa`.
 * Writes results into `out`.
 * Returns 0 on success, -1 if the read cannot be parsed.
 *
 * Thread-safe: stateless.  All working storage is on the stack or in *out.
 */
static int parse_read(const char *seq, const ParseArgs *pa, ParseResult *out)
{
    int slen = (int)strlen(seq);
    int i5_start, i5_end, i3_start, i3_end;

    /* Locate core pattern: [inner5][barcode][inner3] */
    if (find_core(seq, slen, pa, &i5_start, &i5_end, &i3_start, &i3_end) < 0)
        return -1;

    /* ── UMI5: region between outer5 (or read start) and inner5 ── */
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

    /* ── UMI3: region between inner3 and outer3 (or read end) ── */
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

    /* Safety: guard against buffer overflow before copying */
    if (u5_len >= MAX_UMI || pa->barcode_length >= MAX_BC || u3_len >= MAX_UMI)
        return -1;

    memcpy(out->umi5,    seq + u5_start,       (size_t)u5_len);            out->umi5[u5_len] = '\0';
    memcpy(out->barcode, seq + i5_end,  (size_t)pa->barcode_length); out->barcode[pa->barcode_length] = '\0';
    memcpy(out->umi3,    seq + u3_start,       (size_t)u3_len);            out->umi3[u3_len] = '\0';
    return 0;
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

    /* Standard weighted Gini sum, matching the Python implementation */
    double num = 0.0;
    for (int i = 0; i < n; i++)
        num += (double)s[i] * (2.0 * (i + 1) - n - 1);
    free(s);
    return num / ((double)n * (double)total);
}

/* ── FASTQ I/O ──────────────────────────────────────────────────────────── */

/*
 * read_fastq_record: read one 4-line FASTQ record from an open gzFile.
 * All three output buffers must be MAX_LINE bytes.  Trailing whitespace stripped.
 * Returns 1 on success, 0 on EOF or I/O error.
 */
static int read_fastq_record(gzFile fp, char *header, char *seq, char *qual)
{
    char plus[MAX_LINE];
    if (!gzgets(fp, header, MAX_LINE) || header[0] == '\0') return 0;
    if (!gzgets(fp, seq,    MAX_LINE)) return 0;
    if (!gzgets(fp, plus,   MAX_LINE)) return 0;
    if (!gzgets(fp, qual,   MAX_LINE)) return 0;
    header[strcspn(header, "\r\n")] = '\0';
    seq   [strcspn(seq,    "\r\n")] = '\0';
    qual  [strcspn(qual,   "\r\n")] = '\0';
    return 1;
}

/* ── sort comparator ────────────────────────────────────────────────────── */

/* Descending by num_unique_umis (dejackpotted abundance) */
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
        "          [--max-errors <int (default 2)>]\n",
        prog);
    exit(1);
}

/* ── main ───────────────────────────────────────────────────────────────── */

int main(int argc, char **argv)
{
    /* ── argument defaults ── */
    char *input_path  = NULL;
    char *output_path = NULL;
    char *fail_path   = NULL;
    ParseArgs pa = {
        .inner5 = NULL, .inner5_len = 0,
        .inner3 = NULL, .inner3_len = 0,
        .outer5 = NULL, .outer5_len = 0,
        .outer3 = NULL, .outer3_len = 0,
        .barcode_length = 0,
        .min_umi   = 1,
        .max_umi   = 100,
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
            default:  usage(argv[0]);
        }
    }

    if (!input_path || !output_path || !fail_path ||
        !pa.inner5 || !pa.inner3 || pa.barcode_length <= 0)
        usage(argv[0]);

    pa.inner5_len = (int)strlen(pa.inner5);
    pa.inner3_len = (int)strlen(pa.inner3);
    if (pa.outer5) pa.outer5_len = (int)strlen(pa.outer5);
    if (pa.outer3) pa.outer3_len = (int)strlen(pa.outer3);

    /* ── open files ── */
    gzFile in_fp = gzopen(input_path, "r");
    if (!in_fp) { perror(input_path); return 1; }

    FILE *fail_fp = fopen(fail_path, "w");
    if (!fail_fp) { perror(fail_path); return 1; }

    /* ── read loop ──────────────────────────────────────────────────────────
     * NOTE FOR THREADING: this loop body is the natural unit of parallelism.
     * Each worker thread would call parse_read() on its own record and insert
     * into a thread-local `local_tuples` hash table.  After joining all
     * workers, merge the local tables into one global `tuples` table.
     * Writing to fail_fp would require either a per-thread file or a mutex.
     * ──────────────────────────────────────────────────────────────────── */
    TupleEntry *tuples = NULL;
    long n_total = 0, n_parsed = 0, n_failed = 0;
    char header[MAX_LINE], seq[MAX_LINE], qual[MAX_LINE];

    while (read_fastq_record(in_fp, header, seq, qual)) {
        n_total++;
        ParseResult pr;
        int ok = parse_read(seq, &pa, &pr);

        if (ok == 0) {
            char key[MAX_TUPLE_KEY];
            int klen = snprintf(key, sizeof(key), "%s\t%s\t%s",
                                pr.umi5, pr.barcode, pr.umi3);
            /* Treat an overlong key as a failed parse (should not happen with
             * reasonable MAX_UMI / MAX_BC values and checked max_umi). */
            if (klen <= 0 || klen >= (int)sizeof(key)) goto failed;

            TupleEntry *e = NULL;
            HASH_FIND_STR(tuples, key, e);
            if (!e) {
                e = malloc(sizeof(TupleEntry));
                if (!e) { perror("malloc"); return 1; }
                memcpy(e->key, key, (size_t)klen + 1);
                e->count = 0;
                HASH_ADD_STR(tuples, key, e);
            }
            e->count++;
            n_parsed++;
            continue;
        }

    failed:
        n_failed++;
        fprintf(fail_fp, "%s\n%s\n+\n%s\n", header, seq, qual);
    }

    gzclose(in_fp);
    fclose(fail_fp);

    /* ── jackpotting statistics (written to stderr) ── */
    int n_tuples = (int)HASH_COUNT(tuples);
    long *counts_arr = NULL;
    if (n_tuples > 0) {
        counts_arr = malloc((size_t)n_tuples * sizeof(long));
        if (!counts_arr) { perror("malloc"); return 1; }
        int i = 0;
        TupleEntry *e, *tmp;
        HASH_ITER(hh, tuples, e, tmp) counts_arr[i++] = e->count;
    }

    double g = gini(counts_arr, n_tuples);
    long   max_ct  = 0;
    double mean_ct = 0.0;
    for (int i = 0; i < n_tuples; i++) {
        if (counts_arr[i] > max_ct) max_ct = counts_arr[i];
        mean_ct += (double)counts_arr[i];
    }
    if (n_tuples > 0) mean_ct /= n_tuples;
    free(counts_arr);

    fprintf(stderr, "Reads processed:        %ld\n",   n_total);
    fprintf(stderr, "Reads parsed:           %ld (%.1f%%)\n",
            n_parsed, n_total ? 100.0 * n_parsed / n_total : 0.0);
    fprintf(stderr, "Reads failed:           %ld (%.1f%%)\n",
            n_failed, n_total ? 100.0 * n_failed / n_total : 0.0);
    fprintf(stderr, "Unique (UMI,BC) tuples: %d\n",    n_tuples);
    fprintf(stderr, "Gini coefficient:       %.4f  (0 = no jackpotting, 1 = fully jackpotted)\n", g);
    fprintf(stderr, "Mean reads per tuple:   %.1f\n",  mean_ct);
    fprintf(stderr, "Max reads per tuple:    %ld\n",   max_ct);

    /* ── per-barcode summary ── */
    BarcodeSummary *bc_map = NULL;
    {
        TupleEntry *e, *tmp;
        HASH_ITER(hh, tuples, e, tmp) {
            /* Parse key "umi5\tbarcode\tumi3" without strtok */
            const char *tab1 = strchr(e->key, '\t');
            const char *tab2 = tab1 ? strchr(tab1 + 1, '\t') : NULL;
            if (!tab1 || !tab2) continue;

            int bc_len = (int)(tab2 - tab1 - 1);
            char bc_str[MAX_BC];
            if (bc_len <= 0 || bc_len >= MAX_BC) continue;
            memcpy(bc_str, tab1 + 1, (size_t)bc_len);
            bc_str[bc_len] = '\0';

            /* umi-pair key "umi5\tumi3" */
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

            /* add (umi5, umi3) to the per-barcode dedup set */
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

    /* tuple hash no longer needed – free it */
    {
        TupleEntry *e, *tmp;
        HASH_ITER(hh, tuples, e, tmp) { HASH_DEL(tuples, e); free(e); }
    }

    /* umi_set hashes also no longer needed – free them */
    {
        BarcodeSummary *bs, *tmp;
        HASH_ITER(hh, bc_map, bs, tmp) {
            UmiPairEntry *up, *utmp;
            HASH_ITER(hh, bs->umi_set, up, utmp) { HASH_DEL(bs->umi_set, up); free(up); }
        }
    }

    /* collect pointers and sort descending by num_unique_umis */
    int n_bc = (int)HASH_COUNT(bc_map);
    BarcodeSummary **bc_arr = malloc((size_t)n_bc * sizeof(BarcodeSummary *));
    if (!bc_arr && n_bc > 0) { perror("malloc"); return 1; }
    {
        int i = 0;
        BarcodeSummary *bs, *tmp;
        HASH_ITER(hh, bc_map, bs, tmp) bc_arr[i++] = bs;
    }
    qsort(bc_arr, (size_t)n_bc, sizeof(BarcodeSummary *), cmp_bc_desc);

    /* ── write TSV output ── */
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

    /* ── free remaining memory ── */
    {
        BarcodeSummary *bs, *tmp;
        HASH_ITER(hh, bc_map, bs, tmp) { HASH_DEL(bc_map, bs); free(bs); }
    }

    return 0;
}
