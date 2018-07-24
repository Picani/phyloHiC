# -*- coding: utf-8 -*-


"""
genes
=====


This module contains functions that parses BED files containing genes
and compute the adjacency of genes.

A ``gene`` is a ``dict`` with the following keys:

* ``chrom``, a ``str``
* ``start``, an ``int``
* ``end``, an ``int``
* ``name``, a ``str``
* ``strand``, a ``str``, but can only be + or -

In this file, the ``genes`` argument refers to a list of such dicts.

|
:created: May 2018
:last modified: July 2018

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""


def read_bed(name):
    """
    Read the BED file *name*, and return a list of dict. Each dict has the
    following keys: chrom, start, end, name, strand. The list is sorted.

    .. Note:: The BED file is loaded into memory for faster processing.
    """
    with open(name, 'r') as f:
        lines = [l for l in f.read().split('\n') if l]

    res = []
    for line in lines:
        try:
            chrom, start, end, name, _, strand = line.split('\t')
        except ValueError:
            continue
        res.append({
            'name': name,
            'chrom': chrom,
            'start': int(start),
            'end': int(end),
            'strand': strand
        })

    res.sort(key=lambda g: g['start'])
    res.sort(key=lambda g: g['chrom'])
    return res


def make_bed(genes):
    """
    Convert the *genes* to BED. Return the resulting string.
    .. Note:: The `score` column is set to 0.
    """
    rows = [f'{g["chrom"]}\t{g["start"]}\t{g["end"]}\t{g["name"]}\t0\t{g["strand"]}'
            for g in genes]
    return '\n'.join(rows)


def write_bed(genes, name):
    """
    Write the *genes* to the BED file *name*.
    """
    with open(name, 'w') as f:
        f.write(make_bed(genes))
        f.write('\n')


def compute_adjacent(genes):
    """
    Look for the adjacent genes of all genes and return a mapping
    of a gene name to its left and right neighbours. The strand is
    not taken into account for this.

    .. Warning:: *genes* is expected to be sorted.
    """
    res = {}
    # We add the first gene
    res[genes[0]['name']] = {'left': ''}
    if genes[0]['chrom'] == genes[1]['chrom']:
        res[genes[0]['name']]['right'] = genes[1]['name']

    for i, g in enumerate(genes[1:-1]):
        d = {'left': '', 'right': ''}
        if genes[i]['chrom'] == g['chrom']:
            d['left'] = genes[i]['name']
        if genes[i+1]['chrom'] == g['chrom']:
            d['right'] = genes[i+2]['name']

        res[g['name']] = d

    # We add the last gene
    res[genes[-1]['name']] = {'right': ''}
    if genes[-1]['chrom'] == genes[-2]['chrom']:
        res[genes[-1]['name']]['left'] = genes[-2]['name']

    return res
