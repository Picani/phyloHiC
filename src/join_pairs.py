#!/usr/bin/env python3

"""
join_pairs.py
=============

This script reads two pairs files (the left one and the right one)
made with :doc:`make_pairs` and join them based on a list of orthologs.

The result is a data set with 5 columns:

* (str) gene name 1 in species left
* (str) gene name 2 in species left
* (str) gene name 1 in species right
* (str) gene name 2 in species right
* (float) Hi-C value for gene 1 and 2 in species left
* (float) Hi-C value for gene 1 and 2 in species right

This table is written as TSV with no header row, in a Gzipped file.


|
:created: May 2018
:last modified: May 2018

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

import argparse
import logging
import gzip
import csv

from math import isnan
from os.path import splitext
from collections import defaultdict
from distutils.util import strtobool

from iolib import read_orthos


def select_lines(orthos, f):
    """
    Select the records from the pairs file *f* based on the presence of
    their genes in the set *orthos*. *f* shall be already opened.
    Return a dict of that records mapped on their corresponding pair of
    orthology groups.
    """
    reader = csv.reader(f, delimiter='\t')
    # res = {}
    for record in reader:
        if record[0] in orthos and record[1] in orthos:
    #         gr = sorted((groups[record[0]], groups[record[1]]))
    #         res[f'{gr[0]}_{gr[1]}'] = record
    # return res
            yield record


def adjaceny_status(left, right):
    """
    Convert two booleans representing whether the genes are adjacent
    in a string which explicit the relation (either `both` if adjacent
    in both species, `left` if adjacent in left species only, *etc.*)
    """
    if left and right:
        return 'both'
    elif left and not right:
        return 'left'
    elif not left and right:
        return 'right'
    else:
        return 'none'

def threshold(limit, val):
    """
    If *val* < *limit*, return *limit*, else return *val*.
    """
    if val < limit:
        return limit
    else:
        return val

def get_records(records, adj_status, th):
    """
    Get all records from *records* that can be written, filter them based
    on their adjacency, apply the threshold and return a list of tuple
    ready to be written (*i.e.* passed to a csv.writer)

    *records* in a dict of dict of tuple:
    * the first level key is a string of the orthology groups joined
      by a `_`,
    * the second level is either `left` or `right`,
    * the tuple is a record as read by csv.reader from a pairs file.

    The result is a list of tuple corresponding to a row of the file
    described in the above preamble.
    """
    res = []
    ks = list(records.keys())

    done = list()
    for k in ks:
        if len(records[k]) == 2:
            done.append(k)

    for k in done:
        r1 = records[k]['left']
        r2 = records[k]['right']
        del records[k]
        adj = adjaceny_status(bool(strtobool(r1[3])), bool(strtobool(r2[3])))

        if adj_status == "all":
            pass
        elif adj_status == 'none':
            if adj != 'none':
                continue
        elif adj_status == 'and':
            if adj != 'both':
                continue
        elif adj_status == 'or':
            if adj == 'none':
                continue
        elif adj_status == 'xor':
            if adj == 'none' or adj == 'both':
                continue

        if th is None:
            r = (r1[0], r1[1], r2[0], r2[1], r1[2], r2[2])
        else:
            val1 = threshold(th, float(r1[2]))
            if isnan(val1):
                print('dans val1')
                continue
            val2 = threshold(th, float(r2[2]))
            if isnan(val2):
                print('dans val2')
                continue
            r = (r1[0], r1[1], r2[0], r2[1], str(val1), str(val2))

        res.append(r)

    return res



def main():
    desc = ('Reads two pairs files (the left one and the right one) and join '
            'them based on a list of orthologs.')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('orthos', help=('the orthologs file, can be a pair '
                                        'of species or not (see doc)'))
    parser.add_argument('left',
                        help='the left pairs file, optionally Gzipped')
    parser.add_argument('right',
                        help='the right pairs file, optionally Gzipped')
    parser.add_argument('outfile',
                        help='the output file, optionally Gzipped')
    parser.add_argument('-t', '--threshold', type=float,
                        help='the threshold to apply')
    parser.add_argument('-a', '--adjacencies', default='all',
                        choices=['all', 'none', 'and', 'or', 'xor'],
                        help='the adjacencies status to keep')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be verbose')
    args = parser.parse_args()

    if args.verbose:
        level = logging.INFO
    else:
        level = logging.WARN
    logging.basicConfig(format='%(asctime)s: %(message)s', level=level)

    logging.info("C'est parti !")

    orthos, groups = read_orthos(args.orthos)
    logging.info('Loaded orthologs.')

    if splitext(args.left)[1] == '.gz':
        f1 = gzip.open(args.left, 'rt')
    else:
        f1 = open(args.left, 'r')
    lines1 = select_lines(orthos, f1)
    logging.info('Opened left species pairs file.')

    if splitext(args.right)[1] == '.gz':
        f2 = gzip.open(args.right, 'rt')
    else:
        f2 = open(args.right, 'r')
    lines2 = select_lines(orthos, f2)
    logging.info('Opened right species pairs file.')

    if splitext(args.outfile)[1] == '.gz':
        f = gzip.open(args.outfile, 'wt')
    else:
        f = open(args.outfile, 'w')
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')

    records = defaultdict(dict)
    buf = []
    i = 0
    while True:
        try:
            record1 = next(lines1)
            gr = sorted((groups[record1[0]], groups[record1[1]]))
            records[f'{gr[0]}_{gr[1]}']['left'] = tuple(record1)
        except StopIteration:
            record1 = None

        try:
            record2 = next(lines2)
            gr = sorted((groups[record2[0]], groups[record2[1]]))
            records[f'{gr[0]}_{gr[1]}']['right'] = tuple(record2)
        except StopIteration:
            record2 = None

        i += 1
        if i >= 10000:
            i = 0
            buf.extend(get_records(records, args.adjacencies, args.threshold))

        if len(buf) >= 50000:
            writer.writerows(buf)
            buf = []

        if record1 is None and record2 is None:
            break

    logging.info('Done iterating through pairs files.')
    buf.extend(get_records(records, args.adjacencies, args.threshold))
    if buf:
        writer.writerows(buf)

    f1.close()
    f2.close()
    f.close()
    logging.info(f'Written to {args.outfile}')
    logging.info("C'est fini !")


if __name__ == '__main__':
    main()
