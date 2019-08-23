#!/usr/bin/env python3

"""
statshic.py
===========

Compute basic statistics about the given Hi-C dataset.

The statistics are the following:

* Mean
* Standard deviation
* Median
* Percentiles at 10, 25, 75 and 95%

Also report the dimensions of the matrices and their size (number of
elements). For the whole datasets, only the size is given.


:created: August 2018
:last modified: August 2018

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

import json
import argparse
from itertools import combinations_with_replacement

import numpy as np

from hic import HiC, NoSuchHeatmap


def cli_parser():
    desc=('Compute basic statistics about the given Hi-C dataset. '
          'The results are printed to standard output or written as JSON.')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('hic',
                        help='the Hi-C directory with the metadata.json file')
    parser.add_argument('-w', '--whole-dataset', action='store_true',
                        help=('compute also the stats on the whole dataset;'
                              ' may use a lot of memory'))
    parser.add_argument('-o', '--output', help='the output file, in JSON')
    return parser

def main():
    parser = cli_parser()
    args = parser.parse_args()

    hic = HiC(args.hic)

    s = dict()
    if args.whole_dataset:
        whole_dataset = []
    for c1, c2 in combinations_with_replacement(hic.chromosomes, 2):
        try:
            hic.load_map(c1, c2)
        except NoSuchHeatmap:
            continue
        d = {}
        d['dims'] = hic.current['dims']
        d['size'] = d['dims'][0] * d['dims'][1]
        d['mean'] = np.mean(hic.current['data'])
        d['stddev'] = np.std(hic.current['data'])
        # d['stderr'] = d['stddev'] / np.sqrt(len(d))
        d['p10'] = np.percentile(hic.current['data'], 10)
        d['p25'] = np.percentile(hic.current['data'], 25)
        d['median'] = np.median(hic.current['data'])
        d['p75'] = np.percentile(hic.current['data'], 75)
        d['p90'] = np.percentile(hic.current['data'], 90)

        s[f'{c1}|{c2}'] = d
        if args.whole_dataset:
            whole_dataset.extend(hic.current['data'].flatten())

    if args.whole_dataset:
        # w = np.array(whole_dataset)
        w = whole_dataset
        s['all'] = {
            'size': len(w),
            'mean': np.mean(w),
            'stddev': np.std(w),
            # 'stderr'] = 'stddev'] / np.sqrt(len(d))
            'p10': np.percentile(w, 10),
            'p25': np.percentile(w, 25),
            'median': np.median(w),
            'p75': np.percentile(w, 75),
            'p90': np.percentile(w, 90),
        }

    if args.output is None:
        for chroms, d in s.items():
            print(chroms)
            for stat, v in d.items():
                print(f'  {stat} = {v}')
    else:
        with open(args.output, 'w') as f:
            json.dump(s, f)


if __name__ == '__main__':
    main()
