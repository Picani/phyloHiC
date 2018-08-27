#!/usr/bin/env python3

"""
norm_center.py
==============

This script normalizes multiple Hi-C datasets so they can be compared
together. This is done in two steps, each one being done independently
on each dataset.

The first step consists of the following process:

1. All matrices are read into memory.
2. The mean and standard deviation are computed.
3. For each box of the matrices, the mean is subtracted from the box,
   then the box is divided by the standard deviation.
4. The normalized matrices are written.

The second step is simple: we just add to each box in each matrix from
each dataset the minimum value over all dataset.

|
:created: August 2018
:last modified: August 2018

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""


import os
import csv
import sys
import gzip
import logging
import argparse
import concurrent.futures
from shutil import copyfile
from os.path import abspath, basename
from itertools import combinations_with_replacement, product

import numpy as np

from hic import HiC, NoSuchHeatmap


def normalize(hicdata, outdir):
    """
    Normalize the *hicdata* and write them in *outdir*.
    Return the minimum value of the normalized data.
    """
    hicdata.load_all_maps()
    mean = np.mean(hicdata.current['data'])
    stddev = np.std(hicdata.current['data'])
    minimum = float('+inf')

    dd = basename(abspath(outdir))
    logging.debug(f'{dd}: Done loading all maps.')
    try:
        os.mkdir(outdir)
        logging.debug(f'{outdir} created.')
    except FileExistsError:
        pass

    for c1, c2 in combinations_with_replacement(hicdata.chromosomes, r=2):
        try:
            hicdata.load_map(c1, c2)
            logging.debug(f'{dd}: Loaded map for {c1} and {c2}.')
        except NoSuchHeatmap:
            continue

        m = hicdata.current['data'] - mean
        m /= stddev
        minimum = min(minimum, np.min(m))

        rows = []
        for i, j in product(*map(range, hicdata.current['dims'])):
            if m[i, j] == 0.0:
                continue
            rows.append([i*hicdata.binsize, j*hicdata.binsize, m[i, j]])

        c1, c2 = hicdata.current['chroms']
        filename = hicdata._mapfiles[f'{c1}|{c2}']
        logging.debug(f'{dd}: rows prepared for {filename}.')
        with gzip.open(f'{outdir}/{filename}', 'wt') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(rows)
            f.write('\n')
            print(f'{dd}: written {filename}.')

    copyfile(f'{hicdata.datadir}/metadata.json', f'{outdir}/metadata.json')

    logging.info(f'Done normalizing {dd}.')

    return minimum


def add_hic(hicdata, value, verbose=False):
    """
    Substrat *value* from each box of each matrix of the *hicdata*.

    .. Warning:: This is done **in-place**.
    """
    dd = basename(abspath(hicdata.datadir))
    for c1, c2 in combinations_with_replacement(hicdata.chromosomes, r=2):
        try:
            hicdata.load_map(c1, c2)
            logging.debug(f'{dd}: map for {c1} and {c2} loaded.')
        except NoSuchHeatmap:
            continue

        hicdata.current['data'] += value

        rows = []
        for i, j in product(*map(range, hicdata.current['dims'])):
            if hicdata.data['current'][i, j] == 0.0:
                continue
            rows.append([i*hicdata.binsize, j*hicdata.binsize,
                         hicdata.current['data'][i, j]])

        c1, c2 = hicdata.current['chroms']
        filename = hicdata._mapfiles[f'{c1}|{c2}']
        logging.debug(f'{dd}: rows for {filename} prepared.')
        with gzip.open(f'{hicdata.datadir}/{filename}', 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(rows)
            f.write('\n')
            logging.debug(f'{dd}: written {filename}.')

        logging.info(f'Done adding {value} on dataset {dd}')


def cli_parser():
    desc = ('This script normalizes multiple Hi-C datasets so they can be '
            'compared together. This is done by first subtracting the mean '
            'from each box, then by dividing each one by the std. dev. and '
            'finally by adding the last value over all dataset. For more '
            'details, see the head of that file. The normalized datasets are '
            'named after the original ones suffixed with _centernorm')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('datasets', nargs='+', help='the original datasets')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help=('the number of threads to use; '
                              'high values use more memory'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be verbose')
    parser.add_argument('--debug', action='store_true',
                        help='print debug information')

    return parser


def main():
    parser = cli_parser()
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s')
    elif args.verbose:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    else:
        logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(message)s')

    datasets = {}
    norm_datasets = {}
    for d in args.datasets:
        try:
            datasets[d] = HiC(d)
        except Exception as e:
            logging.error(f'Error: {d}: {e}')
            exit(1)
        if d[-1] == '/':
            norm_datasets[d] = f'{d[:-1]}_centernorm'
        else:
            norm_datasets[d] = f'{d}_centernorm'

    logging.info('Normalizing matrices...')

    minimum = float('+inf')
    with concurrent.futures.ProcessPoolExecutor(args.threads) as e:
        futures = []
        for d, h in datasets.items():
            futures.append(e.submit(normalize(h, norm_datasets[d])))

        for future in concurrent.futures.as_completed(futures):
            try:
                data = future.result()
            except Exception as e:
                logging.error(f'Error: {e}')
                exit(1)

            minimum = min(minimum, data)

    logging.info('Subtracting the minimum value...')

    futures = []
    with concurrent.futures.ProcessPoolExecutor(args.threads) as e:
        for d in norm_datasets.values():
            h = HiC(d)
            futures.append(e.submit(add_hic(h, minimum)))

    concurrent.futures.wait(futures)
    logging.info("C'est fini !")


if __name__ == '__main__':
    main()
    sys.exit(0)
