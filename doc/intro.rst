Introduction and General Process
================================

Main Idea
---------

We want to compare the shape of DNA between homologous genomic regions.

The shape of DNA cannot be directly observed on a genome-wide basis, thus we
use the only experiment that offers an approaching information: genome-wide
chromosome conformation capture, also called Hi-C [Lieberman2009]_.

The determination of homologous genomic regions is a complex problem. In order
to simplify it and be able to work more easily with genomes of different size
(around 150Mpb for fly to compare with the multiple Gpb of Human!), we chose to
use orthologous genes. Basically, orthologous genes (or orthologs) are
homologous genes found in different species. For more information on how to get
them, please refer to the `Quest for Orthologs consortium`_ website.

Let *a1* and *b1* be two genes in the genome 1; *a2* and *b2* be two genes in
the genome 2. *a1* and *a2* are orthologs, so are *b1* and *b2*. The idea is
to look at the contacts between *a1* and *b1*, and to compare them with the
contacts between *a2* and *b2*. For that we compute the ratio between their
number of contacts.

To completely compare two genomes, we compute the ratios for all possible pairs
of orthologs.


.. _general_process:

General Process
---------------

The basic process consists of the following steps:

1. **Prepare the data.** Basically, this means selecting the orthologous
   genes that will be used, checking the orthology relationships (only 1-1),
   getting Hi-C contacts, taking care of the Hi-C data normalization.
 
2. **Look at the contacts for all pairs of genes.** This is done inside a
   genome. Thus, this step has to be repeated for each analyzed genome.

3. **Join the contacts between two genomes.** This joining is done by using
   the orthology relationships: the genes *a1* and *b1* in the genome 1 are
   joined to their orthologs counterparts *a2* and *b2* in the genome 2. The
   term *joining* here is just a placeholder for the performed mathematical
   operation ; as said earlier, in our case we just compute a ratio. The
   process is applied to all pairs of orthologs in both genome 1 and 2. This
   step has to be repeated for each pair of analyzed genome.

4. **Make the distances matrix.** This is done by computing a distance
   for each pair of genome using the results of the previous step.

5. **Infer a phylogeny** using the distances matrix computed at the step 4.
   For this step, regular distance-based methods and tools are used.
   Specifically, we used the BioNJ algorithm [Gascuel1997]_ implemented in the
   FastME tool [Lefort2015]_.

Disconnecting steps 2 and 3 enables parallel computation. Since step 2 is
computed on a per-species basis, all species can be run in parallel. Likewise,
all pairs of species can be run in parallel for step 3. in this pipeline, the
parallel computations are handled automagically with Snakemake [Köster2012]_
and will be made explicit in the next sections.

Parameters
----------

There are three main parameters that affect the results.

1. Species and orthologs.
2. Hi-C resolution, also called binsize.
3. Hi-C threshold.

The effect of changing the species is obvious. Changing the orthologs means
looking at different genomic regions. While this will obviously change the
experiment results, we cannot predict in which way. Thus, the whole pipeline
has not been designed to save computation time when changing species and/or
orthologs.

Changing the Hi-C resolution enables more precise results. Indeed, when the
binsize is smaller (and the Hi-C resolution is higher), one can better
distinguish contact locations. This leads to a better map of which gene
contacts which one. However, lowering binsize is a double-edge sword. Since
binning is done in order to get more signal, lowering it dilutes that signal.
This can potentially lead to a loss of information and eventually wrong
results. Finding the good resolution depends on the data and a lot of tries.
In this pipeline, the choice of a binsize is done at step 2. One can run in
parallel multiple steps 2 for the same species, each one with a different
binsize.

The comparison of the contacts for a pair of orthologs between two species is
done using a ratio. We usually don't want to be sensitive to the order of
magnitude of the number of contacts. That is why we chose a ratio. However, if
the contacts from both species are really small, this can mean that they are
not informative (not enough signal). In order to not being affected by this
flaw, a threshold is applied on the contact data during step 3. Testing for
different threshold values can be done by running in parallel different
steps 3 with different threshold values.


.. _Quest for Orthologs consortium: https://questfororthologs.org

.. [Gascuel1997] Gascuel O.,
                 "BIONJ: an improved version of the NJ algorithm based on a
                 simple model of sequence data.",
                 Molecular Biology and Evolution. 1997 14:685-695.

.. [Lieberman2009] Lieberman-Aiden E., van Berkum N. L., Williams L., Imakaev M.,
                   Ragoczy T., Telling A., Amit I., Lajoie B. R., Sabo P. J.,
                   Dorschner M. O., Sandstrom R., Bernstein B., Bender M. A.,
                   Groudine M., Gnirke A., Stamatoyannopoulos J., Mirny L. A.,
                   Lander E. S., Dekker J.
                   "Comprehensive Mapping of Long-Range Interactions Reveals
                   Folding Principles of the Human Genome"
                   Science  09 Oct 2009, Vol. 326, Issue 5950, pp. 289-293
                   doi: 10.1126/science.1181369

.. [Köster2012] Köster J. and Rahmann S.,
                "Snakemake — a scalable bioinformatics workflow engine",
                Bioinformatics, Volume 28, Issue 19, 1 October 2012,
                Pages 2520–2522,
                doi: 10.1093/bioinformatics/bts480

.. [Lefort2015] Lefort V., Desper R. and Gascuel O.,
                "FastME 2.0: A Comprehensive, Accurate, and Fast
                Distance-Based Phylogeny Inference Program",
                Molecular Biology and Evolution, Volume 32, Issue 10,
                1 October 2015, Pages 2798–2800,
                doi: 10.1093/molbev/msv150

