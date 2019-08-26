#!/usr/bin/env python3

"""
dist_all_pairs.py
=================

Compute the scaled distance using all pairs of orthologs (the "standard
way"). The result is **one** distance matrix for a full dataset.

The mode argument define what kind of value we want to keep while computing
the distance:

* `intersection`: only the values present in all species are kept.
* `atLeastTwo`: keep the values present in at least two species.
* `union`: keep all values.

.. note::
   It is intended to be used *instead of* :doc:`/scripts/dist_pairs_indep`
   and :doc:`/scripts/bootstrap`.


:created: August 2019
:last modified: August 2019

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

import argparse
import sys
import re

from toolz import merge_with
from toolz.curried import merge

from distlib import scaled_L2norm, filter_values
from iolib import read_orthos, read_values, phylip


def cli_parser():
    desc = "Compute the distance matrix using all pairs of genes."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('orthos', help='all the orthos, in TSV')
    parser.add_argument('outfile',
                        help='the matrix filename')
    parser.add_argument('values', nargs='+',
                        help='the values files, in (Gzipped) TSV')
    parser.add_argument('-m', '--mode', default='intersection',
                        choices=['intersection', 'atLeastTwo', 'union'],
                        help=('the mode, that is, the kind of values we want'
                              ' to keep while computing the distance'))
    parser.add_argument('-p', '--progress', action='store_true',
                        help='print a progress bar; need tqdm to be installed')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be verbose')
    return parser


def main():
    parser = cli_parser()
    args = parser.parse_args()

    if args.progress:
        try:
            from tqdm import tqdm
            args.verbose = False
        except ImportError:
            print('error: -p/--progress needs tqdm to be installed.')
            sys.exit(1)

    orthos, groups = read_orthos(args.orthos)
    if args.verbose:
        print('Read orthologs.')
    species = set()
    values = dict()

    name_pattern = re.compile(r'(?:\S+/)*(\w+)_(\w+)_values.tsv.gz', re.IGNORECASE)
    if args.progress:
        print('Reading values files...')
        pbar = tqdm(total=len(args.values))
    for filename in args.values:
        m = name_pattern.match(filename)
        if m is None:
            print('Ignored {}'.format(filename))
            continue
        sp_left, sp_right = m.group(1), m.group(2)
        species.add(sp_left)
        species.add(sp_right)
        this_values = read_values(orthos, groups, sp_left, sp_right, filename)
        values = merge_with(merge, values, this_values)
        if args.verbose:
            print(f'Read {filename}')
        elif args.progress:
            pbar.update(1)

    # At this point, we don't need the group number anymore
    values = filter_values(args.mode, species, values.values())

    if args.progress:
        pbar.close()
        print('Computing the distances...')
    elif args.verbose:
        print('Computing the distances...')

    distances = scaled_L2norm(species, values)
    matrix = phylip(species, distances)
    with open(args.outfile, 'w') as f:
        f.write(matrix)
        f.write('\n')


if __name__ == '__main__':
    main()
    sys.exit(0)
