Pairs File
==========

An intermediate file created by the pipeline when looking at the contacts
between each pairs of genes. It is created by the script
:doc:`make_pairs.py </scripts/make_pairs>`.

It is a tabulated-separated values file, usually gzipped. There is no header
row. Each rows corresponds to a pair of genes. Then, the columns are the
following:

* the name for the first gene of the pair, a string,
* the name for the second gene of the pair,  a string,
* the Hi-C value (*i.e.* the corrected and normalized number of contacts) for
  these genes, a floating-point number,
* the adjacency status, *i.e.* if these genes are next to each other along the
  chromosome, a boolean value (either True or False, case-insensitive).
