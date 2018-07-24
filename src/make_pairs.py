#!/usr/bin/env python3

"""
make_pairs.py
=============

Make the pairs of genes for a species.

For each pair, fetch the Hi-C values for two bins where the TSS
are located. Also check the adjacency of the genes. If wanted,
the Hi-C matrices can be scrambled in-memory before use. This
feature uses two functions written by Krister SWENSON for the
**locality** program.

The result is a TSV file, optionally Gzipped (if the file name has
the extension .gz). The file has no header and the following columns:

* gene 1 name (string)
* gene 2 name (string)
* Hi-C value (64 bits float)
* Adjacency status (either True or False, case-insensitive)


|
:created: May 2018
:last modified: July 2018

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

import argparse
import logging
import gzip
import csv
import sys

from math import isnan
from os.path import splitext
from itertools import combinations

import hic
from genes import read_bed, compute_adjacent


def cli_parser():
    parser = argparse.ArgumentParser(
    description='Make the pairs of genes for a species.')
    parser.add_argument('genes', help='the genes, in BED')
    parser.add_argument('hic',
                        help='the directory with the Hi-C sparses matrices')
    parser.add_argument('output', help='explicit enough; can be gzipped')
    parser.add_argument('-N', '--no-nan', action='store_true',
                        help='skip the pair if the Hi-C value is Not-a-Number')
    parser.add_argument('-i', '--intra', action='store_true',
                        help=('only look at genes on the same chromosomes; '
                              'if not set and there is no interchromosomal '
                              'values, Not-a-Number is used'))
    parser.add_argument('-s', '--scramble', action='store_true',
                        help='Scramble in-memory the Hi-C matrices before use')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be verbose')
    return parser


def main():
    parser = cli_parser()
    args = parser.parse_args()

    if args.verbose:
        level = logging.INFO
        # level = logging.DEBUG
    else:
        level = logging.WARN
    logging.basicConfig(format='%(asctime)s: %(message)s', level=level)

    logging.info("C'est parti !")

    genes = read_bed(args.genes)
    adjacencies = compute_adjacent(genes)
    logging.info('Loaded genes.')

    exp = hic.HiC(args.hic)
    logging.info('Loaded Hi-C')

    if splitext(args.output)[1] == '.gz':
        outfile = gzip.open(args.output, 'wt')
    else:
        outfile = open(args.output, 'w')

    writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
    logging.info('Beginning to write pairs...')

    buf = []
    for g1, g2 in combinations(genes, 2):
        c1 = g1['chrom']
        c2 = g2['chrom']

        if c1 not in exp.chromosomes or c2 not in exp.chromosomes:
            continue
        if c1 != c2 and args.intra:
            continue
        if c1 != c2 and not exp.inter:
            value = float('nan')
        else:
            if (c1, c2) != exp.current['chroms'] and\
               (c2, c1) != exp.current['chroms']:
                try:
                    exp.load_map(c1, c2, args.scramble)
                    logging.debug(f'Loaded {c1} x {c2}')
                except hic.NoSuchHeatmap as e:
                    logging.warn(e)
                    continue

            value = exp.get_contact(g1, g2)

        if args.no_nan and isnan(value):
            continue

        a1 = adjacencies[g1['name']]
        adj = a1['left'] == g2['name'] or a1['left'] == g2['name']

        buf.append([g1['name'], g2['name'], str(value), str(adj)])
        if len(buf) == 100000:
            writer.writerows(buf)
            buf = []
            logging.debug('Write')

    if buf:
        writer.writerows(buf)
    outfile.close()
    logging.info("C'est fini !")


if __name__ == '__main__':
    main()
    sys.exit(0)
