#!/usr/bin/env python3


"""
informative_traits.py
=====================

Report information about the traits that are informative, *i.e.* the traits
that don't have the same value in all species.

The result is a key-value association, either displayed in the console or
written to a file as TSV or JSON.

The information reported is:

* The total number of values (keyed ``TotalSize``).
* The number of informative values (keyed ``InformativesSize``).
* The percentage of informative values (keyed ``PercentageInformative``).
* The first non informative value seen during the computation (keyed
  ``NonInformativeValue``); this was used mainly for debugging.


:created: September 2018
:last modified: August 2019

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""


import argparse
import random
import json
import sys
import os
import re

from math import sqrt, isnan
from itertools import combinations
from collections import defaultdict

from toolz import merge_with
from toolz.curried import merge

from distlib import filter_values
from iolib import read_orthos, read_values


def cli_parser():
    desc = "Report information about the traits that are informative"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('orthos', help='all the orthos, in TSV')
    parser.add_argument('values', nargs='+',
                        help='the values files, in (Gzipped) TSV')
    parser.add_argument('-o', '--output',
                        help='write the output there; can be JSON or TSV')
    parser.add_argument('-m', '--mode', default='intersection',
                        choices=['intersection', 'informative', 'union'],
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

    if args.progress:
        pbar.close()

    res = {}

    # At this point, we don't need the group number anymore
    values = filter_values(args.mode, species, values.values())
    res['TotalSize'] = len(values)

    informatives = [v for v in values if len(set(v.values())) != 1]
    res['InformativesSize'] = len(informatives)
    res['PercentageInformative'] = float(len(informatives)) / float(len(values)) * 100.0

    non_informative_value = float('nan')
    for v in values:
        if len(set(v.values())) == 1:
            non_informative_value = next(iter(v.values()))
            break
    res['NonInformativeValue'] = non_informative_value

    if args.output is None:
        for k, v in res.items():
            print(f'{k}: {v}')
    else:
        _, ext = os.path.splitext(args.output)
        if ext == '.json':
            data = json.dumps(res)
        elif ext == '.tsv':
            data = '\n'.join([f'{k}\t{v}' for k, v in res.items()])
        else:
            print(f'Error: unknown output format extension:  {ext}')
            sys.exit(1)

        with open(args.output, 'w') as f:
            f.write(data)
            f.write('\n')


if __name__ == '__main__':
    main()
    sys.exit(0)
