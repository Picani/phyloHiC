Scripts and API documentation
=============================

Scripts
-------

About the pipeline:

1. :doc:`scripts/make_pairs` which computes the contacts of pairs of genes.

2. :doc:`scripts/join_pairs` which extracts the intersection of pairs from two
   tables of contacts. The intersection is computed using the orthologs.

3. From there, we have two possibilities:

   * :doc:`scripts/dist_all_pairs` if we want to look at all pairs,
   * :doc:`scripts/dist_pairs_indep` if we want to look at each pair
     independently.

   
Tools (not exhaustive):

* :doc:`scripts/statshic` computes basic statistics over a dataset and output the
  results directly or write it as JSON.

* :doc:`scripts/norm_center` normalizes a buch of dataset in order to make them
  comparable with each other.

Scripts usage:

.. toctree::
   :maxdepth: 1

   scripts/make_pairs
   scripts/join_pairs
   scripts/bootstrap
   scripts/dist_all_pairs
   scripts/dist_pairs_indep
   scripts/statshic
   scripts/norm_center

API
---

The scripts make use of the following modules:

.. toctree::
   :maxdepth: 1

   api/distlib
   api/genes
   api/hic
   api/iolib


..
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
