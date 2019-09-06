.. PhyloHiC documentation master file, created by
   sphinx-quickstart on Tue Jun  5 11:14:13 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PhyloHiC's documentation!
====================================

PhyloHiC is a set of scripts that allows the computation of distances between
genomes based on *shared genomic regions* and *genomic contact data*.
Currently, we use orthologous genes as shared regions and Hi-C data as contacts
information. Thus, you must already have these data.


License
-------

The source code is available `here`_ and is licensed under the CeCILL V2.1
since it was written when I was working in a `CNRS`_ laboratory. The CeCILL
license is close to the GNU GPL and compatible with it. The full text for
the license is available in the `LICENSE file in the repository`_.


Citing
------

PhyloHiC has been written during my PhD thesis work, and is described in the
third chapter of `my thesis`_ (in French). Moreover, the results presented
there come from the pipeline documented here. Thus, if you use this pipeline,
please cite my thesis:

::

    Sylvain Pulicani.
    Lien entre les réarrangements chromosomiques et la structure de la chromatine chez la Drosophile.
    Autre [cs.OH]. Université Montpellier, 2018. Français.
    NNT: 2018MONTS105
    HAL: tel-02161932

Please note that a scientific publication is coming.


-------------------------------------------------------------------------------


This documentation is divided in three sections:

.. toctree::
   :numbered:
   :maxdepth: 2

   intro
   definitions
   analysis
   scripts_and_api
   file_formats


.. _my thesis: https://tel.archives-ouvertes.fr/tel-02161932
.. _here: https://gite.lirmm.fr/pulicani/phyloHiC
.. _CNRS: https://www.cnrs.fr/en
.. _LICENSE file in the repository: https://gite.lirmm.fr/pulicani/phyloHiC/blob/master/LICENSE
