Values File
===========

An intermediate file created by the pipeline when "joining" the contacts
between different species based on the orthology relationships. It is
created by the script :doc:`join_pairs.py </scripts/join_pairs>`.

It is a tabulated-separated values file, usually gzipped. There is no header
row. For each row, we have two species (named *left* and *right*), and four
genes, two by species. The gene 1 in species left and the gene 1 in species
right are orthologs, as are the gene 2 in species left and the gene 2 in
species right. Then, the columns are the following:

* the name of the gene 1 in species left, a string,
* the name of the gene 2 in species left, a string,
* the name of the gene 1 in species right, a string,
* the name of the gene 2 in species right, a string,
* the Hi-C value (*i.e.* the corrected and normalized number of contacts) for
  the genes 1 and 2 in species left, a floating-point number,
* the Hi-C value for the genes 1 and 2 in species right.

