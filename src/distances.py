# -*- coding: utf-8 -*-


"""
distances
======

This module contains functions to compute phylogenetic distances.

|
:created: July 2018
:last modified: July 2018

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

from collections import defaultdict
from itertools import combinations
from math import isnan, sqrt


class NaNinValuesError(Exception):
    """
    This exception should be thrown by distance-computing functions
    when they encounter a NaN.
    """

    def __init__(self, values):
        super().__init__(f'NaN in the values:\n{values}')


def scaled_ratio_l2(values):
    """
    Compute the distances using the L2-norm of the ratios scaled by
    the number of that ratios. The ratios are made using the *values*.

    Args:
      values(iterable): The pairs of gene values, as an iterable of dicts
                        with species (str) mapped to values (float).

    Return:
      dict of dict: the distances (float) doubly keyed on species.

    .. warning:: This function crashes on Not-a-Number (nan) values.
    """
    distances = defaultdict(dict)
    sizes = defaultdict(dict)

    for v in values:
        for sp1, sp2 in combinations(v.keys(), 2):
            if sp1 not in distances:
                # If sp1 is not in distances, that means:
                # 1. we have never seen sp1 in previous values and we want to
                #    had it in distances.
                # 2. we have seen it but at second-key species, so sp2 is
                #    already in distances.
                sp1, sp2 = sp2, sp1
            if distances[sp1].get(sp2) is None:
                distances[sp1][sp2] = 0.0
                sizes[sp1][sp2] = 0.0

            if isnan(v[sp1]) or isnan(v[sp2]):
                raise NaNinValuesError(v)

            ma = max(v[sp1], v[sp2])
            mi = min(v[sp1], v[sp2])

            if ma == 0.0:
                distances[sp1][sp2] += 0.0  # useless but explicit
            else:
                distances[sp1][sp2] += (1.0 - (mi/ma))**2
            sizes[sp1][sp2] += 1.0

    for sp1 in sizes:
        for sp2 in sizes[sp1]:
            if sizes[sp1][sp2] == 0.0:
                distances[sp1][sp2] = 1.0
            else:
                distances[sp1][sp2] = sqrt(distances[sp1][sp2]) / sqrt(sizes[sp1][sp2])

    return distances


def scaled_count_pairs(values, threshold, above=False):
    """
    Count the number of pairs of genes from *values* for which the *threshold*
    apply. If *above* is True, then the *threshold* applies when the ratio is
    **greater than or equal** to the *threshold*. Else, the *threshold*
    applies when the ratio is **lesser than or equal** to the *threshold*.

    The distances are then scaled to the maximal number of pairs.

    Args:
      values(iterable): The pairs of gene values, as an iterable of dicts
                        with species (str) mapped to values (float).

    Return:
      dict of dict: the distances (float) doubly keyed on species.
    """
    distances = defaultdict(dict)
    sizes = defaultdict(dict)

    for v in values:
        for sp1, sp2 in combinations(v.keys(), 2):
            if sp1 not in distances:
                # Same comment as above.
                sp1, sp2 = sp2, sp1
            if distances[sp1].get(sp2) is None:
                distances[sp1][sp2] = 0.0
                sizes[sp1][sp2] = 0.0

            if isnan(v[sp1]) or isnan(v[sp2]):
                raise NaNinValuesError(v)

            ma = max(v[sp1], v[sp2])
            mi = min(v[sp1], v[sp2])

            if ma == 0.0:
                ratio = 0.0
            else:
                ratio = mi/ma

            if above and ratio >= threshold:
                distances[sp1][sp2] += 1.0
            elif ratio <= threshold:
                distances[sp1][sp2] += 1.0

            sizes[sp1][sp2] += 1.0

    for sp1 in sizes:
        for sp2 in sizes[sp1]:
            if sizes[sp1][sp2] == 0.0:
                distances[sp1][sp2] = 1.0
            else:
                distances[sp1][sp2] = sqrt(distances[sp1][sp2]) / sqrt(sizes[sp1][sp2])

    return distances
