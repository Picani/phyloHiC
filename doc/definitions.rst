A few concepts to define
========================

Before to jump in the next section, let's agree on a few term and their
associated concepts we'll use.

Genes Adjacency
---------------

Two genes are `adjacent` if and only if they are on the same chromosome (or
DNA molecule) and next to each other on the DNA sequence. Obviously, this
concept is only relevant for genes from the same genome.

Now, let :math:`S1` be a set of genes from a species 1 and :math:`S2` be a
set of genes from a species 2. Then let

.. math::
   S1 \times S2 = \{( (a1, b1), (a2, b2) )\}

with

* :math:`a1 \in S1`,
* :math:`b1 \in S1`,
* :math:`a2 \in S2`,
* :math:`b2 \in S2`,
* :math:`a1` and :math:`a2` are orthologous,
* :math:`b1` and :math:`b2` are orthologous.

This corresponds to the step 3 of the :ref:`general_process`.

When collecting the pairs of orthologous from :math:`S1 \times S2`, we
can filter the pairs of pairs of genes based on the adjacency in :math:`S1`
and :math:`S2`. The following situations are possible:

* The genes are adjacent in both pairs (:math:`a1` and :math:`b1` are adjacent
  as are :math:`a2` and :math:`b2`).
* The genes are *not* adjacent in both pairs (:math:`a1` and :math:`b1` are
  *not* adjacent as are :math:`a2` and :math:`b2`).
* The genes are adjacent in :math:`S1` but not in :math:`S2` (:math:`a1` and
  :math:`b1` are adjacent while :math:`a2` and :math:`b2` are *not* adjacent).
* The genes are *not* adjacent in :math:`S1` but are in :math:`S2` (:math:`a1`
  and :math:`b1` are *not* adjacent while :math:`a2` and :math:`b2` are
  adjacent).

We define the following strategies of collecting such pairs of pairs of genes:

* All pairs are collected without taking the adjacency into account; we call
  this case ``all``.
* The pairs are collected only when both genes pairs are *not* adjacent; we
  call this case ``none``.
* The pairs are collected only when the genes of one genes pair are adjacent
  while the genes of the other genes pair are not; we call this case ``xor``.
* The pairs are collected only when at least one of the genes pairs has its
  genes adjacent; we call this case ``or``.
* The pairs are collected only when the genes are adjacent in both genes
  pairs; we call this case ``and``.


Values selection mode
---------------------

The result of the step 3 of the :ref:`general_process` is a list of values
for each pairs of pairs of genes between each pairs of genomes. This list
can contain the same pairs of orthologs among the different pairs of genomes
or not. This depends on the orthologs chosen in the first place and on the
availability of contacts data along each genomes.

The previous facts can have the following effect: different pairs of genomes
raises (partially or totally) different pairs of orthologs. Since we cannot
know *a priori* if this is an issue or not, we put names on the different
situations in order to make it possible to work with each of them in step 4
of the :ref:`general_process`.

These situations are:

* All values are kept, whatever they are present in all or just a subset of
  the pairs of species; we call that situation ``union``.
* Only the values that are present in all pairs of species are kept; we call
  that situation ``intersection``.
* Only the values that are present in at least two pairs of species (so at
  least 3 species) are kept; we soberly call that situation ``atLeastTwo``.
