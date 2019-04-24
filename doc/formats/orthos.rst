Orthologous Genes File
======================

This file contains the orthology relationship between the different species.
It's a tabulated-separated values file, *not* gzipped. The header row
contains the species names or ID (or taxonomic IDs or whatever names you used
for the species in the :doc:`configuration file <config>`). The other rows
contain the gene identifiers used in the genes locations files (those files
are in `BED format`_, so this identifiers are likely the gene IDs from the
database used to get the genes).

On a row, all genes are orthologous together.

Here's an example with six species:
::
 10090              46245       7227        7240        7245        9606
 ENSMUSG00000035948 FBgn0079584 FBgn0039184 FBgn0192545 FBgn0240653 ENSG00000111058
 ENSMUSG00000039632 FBgn0072823 FBgn0036219 FBgn0186022 FBgn0239103 ENSG00000198003
 ENSMUSG00000031578 FBgn0070523 FBgn0030067 FBgn0196092 FBgn0233372 ENSG00000198042
 ENSMUSG00000014856 FBgn0080973 FBgn0034059 FBgn0182903 FBgn0229535 ENSG00000168701
 ENSMUSG00000026869 FBgn0245366 FBgn0030457 FBgn0188676 FBgn0234534 ENSG00000095261
 ENSMUSG00000018845 FBgn0075456 FBgn0010812 FBgn0191449 FBgn0242058 ENSG00000141161
 ENSMUSG00000033629 FBgn0081643 FBgn0032524 FBgn0193470 FBgn0229380 ENSG00000074696
 ENSMUSG00000004264 FBgn0073511 FBgn0010551 FBgn0187108 FBgn0229780 ENSG00000215021
 ENSMUSG00000025939 FBgn0080184 FBgn0033544 FBgn0182517 FBgn0230555 ENSG00000104343


.. _BED format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
