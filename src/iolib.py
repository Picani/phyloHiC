# -*- coding: utf-8 -*-


"""
iolib
======


This module contains functions that parses tabular files used to
store orthologs or pair values, and to format output such as making
PHYLIP matrices.


:created: June 2018
:last modified: August 2019

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

import gzip

from math import isnan


def read_orthos(name):
    """
    Read the orthologs from a TSV file and return a set of the orthologs
    names (from all species) and a mapping of the ortholog names to their
    orthology group.

    The orthologs file is described :doc:`in the documentation </formats/orthos>`.
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

    The values file is described :doc:`in the documentation </formats/values>`.

    .. note::
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


def phylip(species, distances):
    """
    Make the distance matrix for the wanted *species* in PHYLIP
    format. The resulting matrix is returned as a single string
    object.

    .. Warning:: If not all species are present in distances,
                 a KeyError exception is raised.
    """
    lines = [str(len(species))]
    for sp1 in species:
        temp_line = [sp1]
        for sp2 in species:
            if sp1 == sp2:
                temp_line.append('0')
            elif sp1 not in distances or sp2 not in distances[sp1]:
                temp_line.append(f'{distances[sp2][sp1]:.8f}')
            else:
                temp_line.append(f'{distances[sp1][sp2]:.8f}')
        lines.append('\t'.join(temp_line))
    return '\n'.join(lines)
