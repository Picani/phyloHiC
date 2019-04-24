Preparing Data and Running an Analysis
======================================

Here are the steps needed to run one analysis (you can also read one
experiment). By this, we mean measuring the distance between all species
for a fixed set of species, a fixed set of orthologous genes and a fixed set
of "ready-to-use" Hi-C datasets (understand corrected for experimental biases
and normalized). However, the pipeline does allow you to compute the distances
for several Hi-C resolutions/binsizes.

Prepare the data
----------------

1. Prepare the gene locations files for each species, in usual `BED format`_.

2. Prepare the orthologs file in the :doc:`appropriate format <formats/orthos>`.

3. Prepare Hi-C data in the :doc:`appropriate format <formats/hic>`.

4. Prepare the configuration file, in `YAML`_, using the
   :doc:`commented sample <formats/config>`. This file is used by SnakeMake,
   so keep it safely. Using the same config file with the same datasets
   guarantees to re-compute the very same results.

Running the pipeline
--------------------

You just need to launch SnakeMake (with a configuration file named
:file:`config.yaml`):

:command:`snakemake --configfile config.yaml`

|

You can ask SnakeMake to perform multiple steps at once, if possible. For
example, to use 6 jobs at the same time:

:command:`snakemake --configfile config.yaml -j6`

|

SnakeMake can also automagically spans its jobs on a clustering system.
However, be aware that this functionality is system-dependant. Here is a
basic example with an SGE scheduler:

:command:`snakemake --configfile config.yaml -j6 --cluster \'qsub -o outfile -e errfile'`


Results
-------

The pipeline outputs a distance matrix in `PHYLIP format`_ called
:file:`all_replicates.phylip`. It also creates a number of intermediate files
that are kept in case other analysis should be performed. This files are :

* the `pairs files` which are described :doc:`here <formats/pairs>`,
* the `values files` which are described :doc:`here <formats/values>`,
* the `stats file` which is described :doc:`here <formats/stats>`.


.. _BED format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _YAML: https://yaml.org
.. _PHYLIP format: http://evolution.genetics.washington.edu/phylip/doc/distance.html
