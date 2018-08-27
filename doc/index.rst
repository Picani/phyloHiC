.. PhyloHiC documentation master file, created by
   sphinx-quickstart on Tue Jun  5 11:14:13 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PhyloHiC's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


The pipeline
------------

1. :doc:`make_pairs` which computes the contacts of pairs of genes.

2. :doc:`join_pairs` which extracts the intersection of pairs from two tables of contacts. The intersection is computed using the orthologs.

3. From there, we have two possibilities:

   * :doc:`bootstrap` if we want to look at all pairs,
   * :doc:`dist_pairs` if we want to look at each pair independently.

   
Tools
-----

* :doc:`statshic` computes basic statistics over a dataset and output the
  results directly or write it as JSON.


The File Formats
----------------

Coming soon!


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
