#!/usr/bin/env python3


"""
bootstrap.py
============

Randomly sample the pairs of orthologs and compute the scaled distance using
that sample. The sample as the same size as the original set of orthologs
(the sampling is done with replacement). This is the first step of a
bootstrap.

The mode argument define what kind of value we want to keep while computing
the distance:

* `intersection`: only the values present in all species are kept.
* `atLeastTwo`: keep the values present in at least two species.
* `union`: keep all values.

.. note::
   It is intended to be used *instead of* :doc:`/scripts/dist_all_pairs`
   and :doc:`/scripts/dist_pairs_indep`.


:created: May 2018
:last modified: August 2019

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

import argparse
import random
import sys
import os
import re

from toolz import merge_with
from toolz.curried import merge

from distlib import scaled_L2norm, filter_values
from iolib import read_orthos, read_values, phylip


def cli_parser():
    desc = ('Perform the first step of a bootstrap by sampling at random the '
            ' orthologs and computing a distances matrix with that sample.')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('orthos', help='all the orthos, in TSV')
    parser.add_argument('outdir', help='the dir for the resulting samples')
    parser.add_argument('n', type=int, help='the number of replicates')
    parser.add_argument('values', nargs='+',
                        help='the values files, in (Gzipped) TSV')
    parser.add_argument('-o', '--one-file', action='store_true',
                        help='put all matrices in one file called\
                        all_replicates.phylip')
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

    try:
        os.mkdir(args.outdir)
        if args.verbose or args.progress:
            print(f'directory {args.outdir} created.')
    except FileExistsError:
        print('error: {} already exists.'.format(args.outdir))
        sys.exit(1)

    orthos, groups = read_orthos(args.orthos)
    if args.verbose:
        print('Read orthologs.')
    species = set()
    values = dict()

    name_pattern = re.compile('(?:\S+/)*(\w+)_(\w+)_values.tsv.gz', re.IGNORECASE)
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
        print('Computing replicates...')
        pbar = tqdm(total=args.n)

    if args.one_file:
        f = open(args.outdir + '/all_replicates.phylip', 'w')

    for i in range(args.n):
        if args.verbose:
            print(f'Replicate {i}')
        current_choices = random.choices(values, k=len(values))
        distances = scaled_L2norm(species, current_choices)
        matrix = phylip(species, distances)

        if args.one_file:
            f.write(matrix)
            f.write('\n')
        else:
            filename = f'{args.outdir}/replicate_{i}.phylip'
            with open(filename, 'w') as f:
                f.write(matrix)
                f.write('\n')
            if args.verbose:
                print(f'  Written {filename}')

        if args.progress:
            pbar.update(1)

    if args.one_file:
        f.close()
    if args.progress:
        pbar.close()
        print()  # To have a pretty line in the console


if __name__ == '__main__':
    main()
