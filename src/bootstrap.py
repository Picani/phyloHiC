#!/usr/bin/env python3


import argparse
import random
import gzip
import sys
import os
import re

from math import sqrt, isnan
from itertools import combinations
from collections import defaultdict

from toolz import merge_with
from toolz.curried import merge

BUGGED = False

def read_orthos(name):
    """
    Read the orthologs from a TSV file and return a set of the orthologs
    names (from all species) and a mapping of the ortholog names to their
    orthology group.

    The TSV file is expected to have an header row with the species names,
    one orthology group per row, and all species having one and only one
    ortholog per group.
    """
    lines = []
    with open(name, 'r') as f:
        lines = [l for l in f.read().split('\n') if l]
    orthos = set()
    groups = {}
    for i, line in enumerate(lines[1:]):
        fields = line.split('\t')
        for gene in fields:
            orthos.add(gene)
            groups[gene] = i
    return orthos, groups


def read_values(orthos, groups, sp_left, sp_right, name):
    """
    Read the value file at *name* and return a dict of dict of group ids
    to species name to Hi-C value.

    The value file is the one output by phyloHiC compare, i.e. a
    Gzipped TSV file, with no header row and 6 columns:

    * Species left, gene for ortholog group 1
    * Species left, gene for ortholog group 2
    * Species right, gene for ortholog group 1
    * Species right, gene for ortholog group 2
    * Hi-C value for the pair gene 1/gene 2 in Specie left
    * Hi-C value for the pair gene 1/gene 2 in Specie right

    Extract only the values where all four genes are present in *orthos*.
    """
    lines = []
    with gzip.open(name, 'rt') as f:
        lines = [l for l in f.read().split('\n') if l]

    values = {}
    for line in lines[1:]:
        g1l, g1r, g2l, g2r, v1, v2 = line.split('\t')
        if g1l not in orthos or\
           g1r not in orthos or\
           g2l not in orthos or\
           g2r not in orthos:
            continue

        group1 = set([groups[g1l], groups[g1r]])
        group2 = set([groups[g2l], groups[g2r]])
        if group1 != group2:
            print('AAAAAh! {} {} {}\t{} {} {}'.format(g1l, g1r, group1,
                                                      g2l, g2r, group2))
            continue

        f1, f2 = float(v1), float(v2)
        if isnan(f1) or isnan(f2):
            continue

        gr = sorted([groups[g1l], groups[g1r]])
        values[f'{gr[0]}_{gr[1]}'] = {sp_left: f1, sp_right: f2}

    return values


def distances_scaledl2(species, values):
    """
    Compute the distances using the L2-norm scaled by the size with the
    *values* for all pairs of *species*.

    .. note:: This function makes no assumption of the way we want to
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
            if ma == 0.0 and BUGGED:
                # To mimic the results of phyloHiC compare from commit
                # f0415951e58318460c590aa9a818113f32ad37f0
                distances[sp1][sp2] += 1.0
            elif ma == 0.0:
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


def phylip(species, distances):
    lines = [str(len(species))]
    for sp1 in species:
        temp_line = [sp1]
        for sp2 in species:
            if sp1 == sp2:
                temp_line.append(0)
            elif sp1 not in distances or sp2 not in distances[sp1]:
                temp_line.append(distances[sp2][sp1])
            else:
                temp_line.append(distances[sp1][sp2])
        lines.append('\t'.join(map(str, temp_line)))
    return '\n'.join(lines)


def main():
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
    parser.add_argument('-r', '--reference', action='store_true',
                        help=('compute the actual distances matrix and write'
                              ' it in the outdir as ref.phylip'))
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
    args = parser.parse_args()

    global BUGGED
    BUGGED = args.bugged

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
    for i in range(args.n):
        if args.verbose:
            print(f'Replicate {i}')
        current_choices = random.choices(values, k=len(values))
        distances = distances_scaledl2(species, current_choices)

        filename = f'{args.outdir}/replicate_{i}.phylip'
        with open(filename, 'w') as f:
            f.write(phylip(species, distances))
            f.write('\n')
        if args.verbose:
            print(f'  Written {filename}')
        elif args.progress:
            pbar.update(1)

if __name__ == '__main__':
    main()

