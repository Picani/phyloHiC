Hi-C dataset format
===================

An Hi-C dataset is a folder containing the heatmaps and a
:file:`metadata.json` file.


Metadata file
-------------

As its name suggests, the JSON file holds basic metadata about the dataset.
Its format is the following::

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

.. note::
   The `MapFiles` maps to the key-value store of the heatmaps files. The files
   can be named in any way, but not the key. They must be of the form `X|Y`.

Heatmap files
-------------

There are one file per matrix. The matrices are written as tabulated-separated
files, optionally gzipped. There is no header row, and 3 columns:

* row position (an integer),
* column position (an integer),
* value (a floating-point number).

In order to save space, only the non-empty boxes from the matrices are written
(the matrices are sparse).
