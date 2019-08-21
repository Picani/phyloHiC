# -*- coding: utf-8 -*-


"""
distlib
=========

This module contains functions that compute distances and helpers.

:created: August 2019
:last modified: August 2019

.. codeauthor::
   Sylvain PULICANI <pulicani@lirmm.fr>
"""

from collections import defaultdict
from itertools import combinations
from math import isnan, sqrt


def scaled_L2norm(species, values):
    """
    Compute the distances using the L2-norm scaled by the size with the
    *values* for all pairs of *species*.
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
                raise ValueError(f"There shouldn't be any NaN - sp1: {sp1}, sp2: {sp2}")
            ma = max(v[sp1], v[sp2])
            mi = min(v[sp1], v[sp2])
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


def filter_values(constraint, species, values):
    """
    Filter the *values* based upon the *species* and the wanted *constraint*:

    * `intersection`: only the values present in all species are kept.
    * `informative`: keep the values present in at least two species.
    * `union`: keep all values.

    .. note::
       This function is a generator.

    .. warning::
       If *constraint* is something else than `intersection`, `informative`
       or `union`, then a ValueError exception is raised.
    """
    if constraint not in ['intersection', 'informative', 'union']:
        raise ValueError("constraint must be one of 'intersection', 'informative' or 'union'.")

    for value in values:
        if constraint == 'intersection':
            if len(value) == len(species):
                yield value
        elif constraint == 'informative':
            if len(value) != 1:
                yield value
        else:  # in 'union' we keep all values
            yield value

