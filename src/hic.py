# hic.py

import random
import gzip
import json
import csv

from itertools import product, combinations_with_replacement

import numpy as np

class NoSuchHeatmap(Exception):
    pass


class HiC:
    """
    HiC represents an Hi-C experiment, handling the matrices reading,
    and position fetching.

    The experiment is encapsulated into a directory, with one file per matrix.
    The matrices are written as **gzipped TSV** files with 3 columns: `row`
    position, `column` position and `value`. The matrices are **sparse**.

    There is also a JSON file called `metadata.json` of the following form:
    ::
        {
          "Binsize": 5000,
          "Assembly": "droYak2",
          "Species": "Drosophila yakuba",
          "Comment": "nm_none - No NaN",
          "Date": "",
          "Dataset": "Dyak_c",
          "Dims": {
            "2R|2R": [1234, 1234],
            "4|X": [345, 2345],
            ...
          }
          "MapFiles": {
            "2R|2R": "2R_2R.tsv.gz",
            "4|X": "4_X.tsv.gz",
            ...
          }
        }

    Please note that the format of the file names is free, but the one of
    the keys is not. This is of the form `chromosome|chromosome`.

    |
    :created: May 2018
    :last modified: August 2018

    .. codeauthor::
       Sylvain PULICANI <pulicani@lirmm.fr>
    """

    def __init__(self, dirname):
        self.datadir = dirname

        try:
            with open(f'{dirname}/metadata.json', 'r') as f:
                metadata = json.load(f)
        except FileNotFoundError:
            raise Exception('the Hi-C directory does not contain a "metadata.json" file')

        self.binsize = metadata['Binsize']
        self._dims = metadata['Dims']
        self._mapfiles = metadata['MapFiles']
        self.inter = False

        temp = set()
        for chroms in self._mapfiles:
            c1, c2 = chroms.split('|')
            if c1 != c2:
                self.inter = True
            temp.add(c1)
            temp.add(c2)
        self.chromosomes = list(temp)

        self.current = {'data': None, 'dims': (0, 0), 'chroms': ('', '')}


    def load_map(self, rowChrom, colChrom, scramble=False):
        """
        Load the wanted matrix.
        If needed, the row and column chromosomes will be swapped.
        If there is no matrix with this pair of chromosomes, a
        `NoSuchHeatmap` error is raised.
        if *scramble* is `True`, then the matrix is scramble in-place
        after being loaded.
        """
        k = f'{rowChrom}|{colChrom}'
        self.current['chroms'] = (rowChrom, colChrom)
        if k not in self._mapfiles:
            k = f'{colChrom}|{rowChrom}'
            self.current['chroms'] = (colChrom, rowChrom)
        if k not in self._mapfiles:
            raise NoSuchHeatmap(f'{rowChrom} x {colChrom}')

        with gzip.open(f'{self.datadir}/{self._mapfiles[k]}', 'rt') as f:
            data = list(csv.reader(f, delimiter='\t'))

        # max_row_pos = max(set(int(d[0]) for d in data))
        # max_col_pos = max(set(int(d[1]) for d in data))

	    # # Here, the +1 is because we can have something like the last
        # # position read is 11 and the binsize is 5, then we want to
        # # have the following bins:
        # # bin 0 [0, 5), bin 1 [5, 10), bin 2 [10, 15)
	    # # There are 11 / 5 + 1 = 3 bins.
        # nrow = (max_row_pos // self.binsize) + 1
        # ncol = (max_col_pos // self.binsize) + 1
        nrow, ncol = self._dims[k]
        self.current['dims'] = (nrow, ncol)

        m = np.zeros((nrow, ncol), dtype=float)
        for rpos, cpos, value in data:
            ri = int(rpos) // self.binsize
            ci = int(cpos) // self.binsize
            m[ri, ci] = float(value)

        self.current['data'] = m
        if scramble:
            if rowChrom == colChrom:
                self._scrambleIntraFY()
            else:
                self._scrambleInterFY()


    def load_all_maps(self):
        """
        Load all the matrices. The result is an array with all the values.

        .. warning:: This function can need a lot of memory.

        .. warning:: This function's result may change in the future.
        """
        data = []
        for k in self._mapfiles.keys():
            c1, c2 = k.split('|')
            self.load_map(c1, c2)
            data.extend(self.current['data'].flatten())
        self.current['data'] = np.array(data)
        self.current['dims'] = self.current['data'].shape
        self.current['chroms'] = ('all', 'all')


    def get_contact(self, g1, g2):
        """
        Fetch the contact between the genes *g1* and *g2*.
        They are dict of the same kind as returned in the :doc:`genes`.
        If the genes are not on the same chromosomes as the currently loaded
        heatmap, an exception is raised. If at least one of them is outbound,
        Not-a-Number (NaN) is returned.
        """
        if self.current['data'] is None:
            raise Exception('No heatmap is currently loaded')

        c1 = g1['chrom']
        c2 = g2['chrom']
        p1 = g1['start'] // self.binsize
        p2 = g2['start'] // self.binsize

        if (c1, c2) != self.current['chroms']:
            c1, c2 = c2, c1
            p1, p2 = p2, p1

        if (c1, c2) != self.current['chroms']:
            raise Exception((f'{c1} x {c2}: heatmap not loaded; current is '
                             f'{self.current["chroms"][0]} x '
                             f'{self.current["chroms"][1]}'))

        if p1 >= self.current['dims'][0] or p2 >= self.current['dims'][1]:
            # TODO: shouldn't this be NaN?
            # return -1.0
            return float('nan')

        return self.current['data'][p1, p2]


    def _scrambleInterFY(self):
        """
        Randomize the matrix by the 'Fisher-Yates shuffle '<https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#The_modern_algorithm>.

        .. codeauthor::
           Krister SWENSON <swenson@lirmm.fr>
        """
        mat = self.current['data']
        rows = mat.shape[0]
        cols = mat.shape[1]
        for x,y in product(range(rows), range(cols)):
            if mat[x][y] is not None:
                i = random.randrange(rows)
                j = random.randrange(cols)
                while mat[i][j] is None:
                    i = random.randrange(rows)
                    j = random.randrange(cols)

                mat[x][y], mat[i][j] = mat[i][j], mat[x][y]


    def _scrambleIntraFY(self):
        """
        Randomize the matrix by the `Fisher-Yates shuffle `<https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#The_modern_algorithm>.

        .. codeauthor::
           Krister SWENSON <swenson@lirmm.fr>
        """
        mat = self.current['data']
        rows = mat.shape[0]
        cols = mat.shape[1]
        assert(rows == cols)

        for x,y in combinations_with_replacement(range(rows), 2):
            if mat[x][y] is not None:
                #Test for symmetry.
                if mat[x][y] != mat[y][x] and abs(mat[x][y]-mat[y][x]) > .001:
                    raise intraSymmetry(x,y)

                #Find value to swap with.
                i = random.randrange(rows)
                j = random.randrange(rows)
                while(mat[i][j] is None):
                    i = random.randrange(rows)
                    j = random.randrange(cols)

                #Do the swap.
                vals = mat[x][y], mat[i][j]
                mat[i][j], mat[x][y] = vals
                mat[j][i], mat[y][x] = vals
