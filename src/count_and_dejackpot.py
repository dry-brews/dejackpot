'''
This script is designed to count the occurences of barcodes in a sequencing dataset.
The reads are presumed to have the following structure:
[optional constant sequence] - [5' pseudoUMI] - [constant sequence] - [barcode] - [constant sequence] - [3' pseudoUMI] - [optional constant sequence]

Typically, when the combined read lengths is shorter than the amplicon size, the optional constant sequences will not be present after flashing.
When the combined read lengths are longer, flash must be run in --allow-outies mode, and the optional constant sequences will be present after flashing.

The key user-configurable parameters are:
--outer5: the 5' optional constant sequence
--inner5: the 5' constant sequence
--inner3: the 3' constant sequence
--outer3: the 3' optional constant sequence
--barcode-length: the length of the barcode sequence
--max-errors: the maximum allowed number of substitutions per constant sequence match (default: 2)
--min-umi: minimum pseudoUMI length in bp
--max-umi: maximum pseudoUMI length in bp

Usage:
    python count_and_dejackpot.py --in reads.fastq.gz --out counts.tsv --fail unparseable.fastq \
        --inner5 CTGTCTCTTATACA --inner3 ATGTGTGAGAATAG \
        --barcode-length 20 --min-umi 5 --max-umi 15

Output (TSV, stdout or --out):
    barcode  total_reads  num_unique_umis

Jackpotting summary is written to stderr.
'''

import argparse
import gzip
import sys
from collections import defaultdict

import regex


# ---------------------------------------------------------------------------
# FASTQ reading
# ---------------------------------------------------------------------------

def fastq_reader(filepath):
    '''Yield (header, sequence, quality) for each read. Handles .gz files.'''
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq  = f.readline().strip()
            f.readline()            # discard '+' line
            qual = f.readline().strip()
            yield header.strip(), seq, qual


# ---------------------------------------------------------------------------
# Read parsing
# ---------------------------------------------------------------------------

def parse_read(read, core_pattern, outer5_pattern, outer3_pattern, min_umi, max_umi):
    '''
    Extract (umi5, barcode, umi3) from a single read sequence.

    The core_pattern matches: [inner5][barcode][inner3]
    UMI5 is extracted from the region before inner5 (bounded by outer5 if given).
    UMI3 is extracted from the region after inner3 (bounded by outer3 if given).

    Returns None if the read cannot be parsed or UMI lengths are out of bounds.
    '''
    m = core_pattern.search(read)
    if m is None:
        return None

    inner5_start = m.start('inner5')
    inner3_end   = m.end('inner3')
    barcode      = m.group('barcode')

    # UMI5: sequence between outer5 (or read start) and inner5
    if outer5_pattern:
        o5 = outer5_pattern.search(read, 0, inner5_start)
        if o5 is None:
            return None
        umi5 = read[o5.end() : inner5_start]
    else:
        umi5 = read[:inner5_start]

    # UMI3: sequence between inner3 and outer3 (or read end)
    if outer3_pattern:
        o3 = outer3_pattern.search(read, inner3_end)
        if o3 is None:
            return None
        umi3 = read[inner3_end : o3.start()]
    else:
        umi3 = read[inner3_end:]

    # Validate UMI lengths
    if not (min_umi <= len(umi5) <= max_umi):
        return None
    if not (min_umi <= len(umi3) <= max_umi):
        return None

    return umi5, barcode, umi3


# ---------------------------------------------------------------------------
# Jackpotting assessment
# ---------------------------------------------------------------------------

def gini(values):
    '''
    Gini coefficient over a list of non-negative counts.
    Returns 0 for a perfectly even distribution, approaching 1 for extreme inequality.
    '''
    n = len(values)
    if n == 0:
        return float('nan')
    s = sorted(values)
    total = sum(s)
    if total == 0:
        return 0.0
    cumsum = sum(x * (2 * (i + 1) - n - 1) for i, x in enumerate(s))
    return cumsum / (n * total)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Count barcodes in a FASTQ file using pseudoUMIs to detect PCR jackpotting.'
    )
    parser.add_argument('--in',   dest='input',   required=True,
                        help='Input FASTQ or FASTQ.gz file')
    parser.add_argument('--out',  dest='output',  required=True,
                        help='Output TSV file (barcode counts)')
    parser.add_argument('--fail', dest='fail',    required=True,
                        help='Output FASTQ file for reads that could not be parsed')
    parser.add_argument('--inner5', required=True,
                        help="Inner 5' constant sequence (flanks pseudoUMI5 on its 3' side)")
    parser.add_argument('--inner3', required=True,
                        help="Inner 3' constant sequence (flanks pseudoUMI3 on its 5' side)")
    parser.add_argument('--outer5', default=None,
                        help="Outer 5' constant sequence (optional; defines where pseudoUMI5 begins)")
    parser.add_argument('--outer3', default=None,
                        help="Outer 3' constant sequence (optional; defines where pseudoUMI3 ends)")
    parser.add_argument('--barcode-length', dest='barcode_length', type=int, required=True,
                        help='Barcode length in bp')
    parser.add_argument('--min-umi', dest='min_umi', type=int, default=1,
                        help='Minimum pseudoUMI length in bp')
    parser.add_argument('--max-umi', dest='max_umi', type=int, default=100,
                        help='Maximum pseudoUMI length in bp')
    parser.add_argument('--max-errors', dest='max_errors', type=int, default=2,
                        help='Maximum substitutions allowed per constant sequence match (default: 2)')
    args = parser.parse_args()

    # Compile patterns once before the read loop.
    # BESTMATCH selects the alignment with the fewest total substitutions.
    core_pattern = regex.compile(
        rf'(?P<inner5>{args.inner5}){{e<={args.max_errors}}}'
        rf'(?P<barcode>.{{{args.barcode_length}}})'
        rf'(?P<inner3>{args.inner3}){{e<={args.max_errors}}}',
        flags=regex.BESTMATCH
    )
    outer5_pattern = (
        regex.compile(rf'(?:{args.outer5}){{e<={args.max_errors}}}')
        if args.outer5 else None
    )
    outer3_pattern = (
        regex.compile(rf'(?:{args.outer3}){{e<={args.max_errors}}}')
        if args.outer3 else None
    )

    # counts maps (umi5, barcode, umi3) -> number of reads with that exact combination
    counts = defaultdict(int)
    n_total = n_parsed = n_failed = 0

    with open(args.fail, 'w') as fail_fq:
        for header, seq, qual in fastq_reader(args.input):
            n_total += 1
            result = parse_read(seq, core_pattern, outer5_pattern, outer3_pattern,
                                args.min_umi, args.max_umi)
            if result is not None:
                counts[result] += 1
                n_parsed += 1
            else:
                n_failed += 1
                fail_fq.write(f'{header}\n{seq}\n+\n{qual}\n')

    # --- Jackpotting assessment (stderr) ------------------------------------
    all_counts = list(counts.values())
    g = gini(all_counts)
    max_count = max(all_counts, default=0)
    mean_count = sum(all_counts) / len(all_counts) if all_counts else 0

    print(f'Reads processed:        {n_total:,}',            file=sys.stderr)
    print(f'Reads parsed:           {n_parsed:,} ({100*n_parsed/max(n_total,1):.1f}%)', file=sys.stderr)
    print(f'Reads failed:           {n_failed:,} ({100*n_failed/max(n_total,1):.1f}%)', file=sys.stderr)
    print(f'Unique (UMI,BC) tuples: {len(all_counts):,}',   file=sys.stderr)
    print(f'Gini coefficient:       {g:.4f}  (0 = no jackpotting, 1 = fully jackpotted)', file=sys.stderr)
    print(f'Mean reads per tuple:   {mean_count:.1f}',       file=sys.stderr)
    print(f'Max reads per tuple:    {max_count:,}',          file=sys.stderr)

    # --- Per-barcode summary -----------------------------------------------
    # For each barcode, count total reads (before dejackpotting) and the
    # number of unique (umi5, umi3) combinations (the dejackpotted count).
    bc_total_reads  = defaultdict(int)
    bc_unique_umis  = defaultdict(set)

    for (umi5, barcode, umi3), count in counts.items():
        bc_total_reads[barcode] += count
        bc_unique_umis[barcode].add((umi5, umi3))

    # Sort by unique UMI count descending (dejackpotted abundance)
    barcodes_sorted = sorted(
        bc_total_reads.keys(),
        key=lambda bc: len(bc_unique_umis[bc]),
        reverse=True
    )

    with open(args.output, 'w') as tsv:
        tsv.write('barcode\ttotal_reads\tnum_unique_umis\n')
        for bc in barcodes_sorted:
            tsv.write(f'{bc}\t{bc_total_reads[bc]}\t{len(bc_unique_umis[bc])}\n')

    print(f'Unique barcodes:        {len(bc_total_reads):,}', file=sys.stderr)
    print(f'Results written to:     {args.output}',           file=sys.stderr)
    print(f'Failed reads written to: {args.fail}',            file=sys.stderr)


if __name__ == '__main__':
    main()
