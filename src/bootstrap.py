#!/usr/bin/env python3


"""
bootstrap.py
============

TODO: Write this...

"""


import argparse
import random
import sys
import os
import re

from math import sqrt, isnan
from itertools import combinations
from collections import defaultdict

from toolz import merge_with
from toolz.curried import merge

from iolib import read_orthos, read_values, phylip

# BUGGED = False


def distances_scaledl2(species, values):
    """
    Compute the distances using the L2-norm scaled by the size with the
    *values* for all pairs of *species*.

    .. note:: This function makes no assumption on the way we want to
              compute the bootstrap (union or intersection). The caller
              is responsible for that.
    """
    distances = defaultdict(dict)
    sizes = defaultdict(dict)
    species_pairs = list(combinations(species, 2))
    for sp1, sp2 in species_pairs:
        distances[sp1][sp2] = 0.0
        sizes[sp1][sp2] = 0.0

    for v in values:
        for sp1, sp2 in species_pairs:
            if sp1 not in v or sp2 not in v:
                continue
            if isnan(v[sp1]) or isnan(v[sp2]):
                print(v)
                sys.exit(7)
            ma = max(v[sp1], v[sp2])
            mi = min(v[sp1], v[sp2])
            # if ma == 0.0 and BUGGED:
            #     # To mimic the results of phyloHiC compare from commit
            #     # f0415951e58318460c590aa9a818113f32ad37f0
            #     distances[sp1][sp2] += 1.0
            # elif ma == 0.0:
            if ma == 0.0:
                distances[sp1][sp2] += 0.0  # useless but explicit
            else:
                distances[sp1][sp2] += (1.0 - (mi/ma))**2
            sizes[sp1][sp2] += 1.0

    for sp1, sp2 in species_pairs:
        if sizes[sp1][sp2] == 0.0:
            distances[sp1][sp2] = 1.0
        else:
            distances[sp1][sp2] = sqrt(distances[sp1][sp2]) / sqrt(sizes[sp1][sp2])
    return distances


def cli_parser():
    desc = ('Perform the bootstrap (and optionally recompute the actual '
            'distances matrix) from the output values of the "compare" '
            'command of phyloHiC. WARNING: for now, it silently ignores '
            'Not a Number (nan) values.')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('orthos', help='all the orthos, in TSV')
    parser.add_argument('outdir', help='the dir for the resulting samples')
    parser.add_argument('n', type=int, help='the number of replicates')
    parser.add_argument('values', nargs='+',
                        help='the values files, in (Gzipped) TSV')
    parser.add_argument('-o', '--one-file', action='store_true',
                        help='put all matrices in one file called\
                        all_replicates.phylip')
    parser.add_argument('-r', '--reference', action='store_true',
                        help=('compute the actual distances matrix and write'
                              ' it in the outdir as ref.phylip'))
    parser.add_argument('-i', '--informative', action="store_true",
                        help='keep only informative pairs')
    parser.add_argument('-u', '--union', action='store_true',
                        help=('by default, the computation is done with '
                              'the intersection of the orthologs, i.e. if '
                              'a pair of orthologs has no value in at least '
                              'one species, this pair of orthologs is '
                              'discarded for all species; this option removes '
                              'that constraint.'))
    parser.add_argument('-b', '--bugged', action='store_true',
                        help=('a bug was present at the time I generated the '
                              'data I usually use now; this option reproduces '
                              'this bug by putting a 1 instead of a 0 for the '
                              'the pairs of orthologs having both an Hi-C '
                              'value of 0'))
    parser.add_argument('-p', '--progress', action='store_true',
                        help='print a progress bar; need tqdm to be installed')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be verbose')
    return parser


def main():
    parser = cli_parser()
    args = parser.parse_args()

    # global BUGGED
    # BUGGED = args.bugged

    if args.progress:
        try:
            from tqdm import tqdm
            args.verbose = False
        except ImportError:
            print('error: -p/--progress needs tqdm to be installed.')
            sys.exit(1)

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
    values = list(values.values())

    if not args.union:
        values = [v for v in values if len(v) == len(species)]

    if args.informative:
        values = [v for v in values if len(set(v.values)) != 1]

    if args.reference:
        distances = distances_scaledl2(species, values)
        filename = f'{args.outdir}/ref.phylip'
        with open(filename, 'w') as f:
            f.write(phylip(species, distances))
            f.write('\n')
        if args.verbose:
            print(f'Reference matrix written in {filename}')

    if args.progress:
        print('Computing replicates...')
        pbar = tqdm(total=args.n)

    if args.one_file:
        f = open(args.outdir + '/all_replicates.phylip', 'w')

    for i in range(args.n):
        if args.verbose:
            print(f'Replicate {i}')
        current_choices = random.choices(values, k=len(values))
        distances = distances_scaledl2(species, current_choices)
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
        print()  # To have a pretty line in the console


if __name__ == '__main__':
    main()

