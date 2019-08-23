#!/usr/bin/env python3

"""
dist_pairs_indep.py
===================

This script computes one distance matrix by pair of orthologs.
Due to the fact that we want one matrix per pair of genes, this
script runs in intersection mode (*i.e.* only the pairs present
in all species are used).

.. note::
   It is intended to be run *instead of* :doc:`/scripts/dist_all_pairs`
   and :doc:`/scripts/bootstrap`.


:created: June 2018
:last modified: August 2019

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""


import argparse
import sys
import os
import re

from toolz import merge_with
from toolz.curried import merge

from distlib import scaled_L2norm, filter_values
from iolib import read_orthos, read_values, phylip


def cli_parser():
    desc = "Compute the distance matrix on each pair of genes independently."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('orthos', help='all the orthos, in TSV')
    parser.add_argument('outdir',
                        help='the directory in which put the results')
    parser.add_argument('values', nargs='+',
                        help='the values files, in (Gzipped) TSV')
    parser.add_argument('-o', '--one-file', action='store_true',
                        help='put all matrices in one file; outdir is this\
                        file name')
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

    if not args.one_file:
        try:
            os.mkdir(args.outdir)
            if args.verbose:
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
    filter_values('intersection', species, values.values())

    if args.progress:
        pbar.close()
        print('Computing the distances...')
        pbar = tqdm(total=len(values))
    elif args.verbose:
        print('Computing the distances...')

    if args.one_file:
        f = open(args.outdir, 'w')

    for i, v in enumerate(values):
        distances = scaled_L2norm(species, [v])
        matrix = phylip(species, distances)

        if args.one_file:
            f.write(matrix)
            f.write('\n')
        else:
            filename = f'{args.outdir}/matrix_{i}.phylip'
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
    sys.exit(0)
