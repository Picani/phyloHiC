phyloHiC
========

Make a phylogeny using the chromatin structure information extracted from
Hi-C experiments.


Requirements
------------

The code is mainly written using the [Python][8] programming language, with
(currently) one script written using the [Go][9] language and the [gonum][10]
library.

In order to run the pipeline, you will need:

* Python 3 (at least 3.6 for the [f-string][1] notation)
* [PyToolz][2]
* [numpy][3]
* [tqdm][4] (optional, used to display progress bar)


Documentation
-------------

The documentation is available online at [phylohic.rtfd.io][11].

The documentation is written using [sphinx][5]. If you want to build it,
you will need sphinx of course, with the [Read the Doc theme][6] and the
[sphinx-argparse][7] extension.


License
-------

This code has been written during my work as a PhD student in a [CNRS][12]
laboratory. Thus, it is licensed under the CeCILL V2.1, which very close and
compatible with the GNU GPL.

The English text is in the file `LICENSE` and the French text is in the file
`LICENSE_fr`.


Citing
------

This code implements a method described in the third chapter of
[my PhD thesis][13]. Moreover, the results presented there are obtained using
the code in this repository.

If you use it, please cite my thesis:

    Sylvain Pulicani.
    Lien entre les réarrangements chromosomiques et la structure de la chromatine chez la Drosophile.
    Autre [cs.OH]. Université Montpellier, 2018. Français.
    NNT: 2018MONTS105
    HAL: tel-02161932

Please note that a scientific publication is coming.


[1]: https://www.python.org/dev/peps/pep-0498/
[2]: https://toolz.readthedocs.io/en/latest/
[3]: http://www.numpy.org/
[4]: https://pypi.org/project/tqdm/
[5]: http://www.sphinx-doc.org/en/master/
[6]: https://github.com/rtfd/sphinx_rtd_theme
[7]: https://github.com/ribozz/sphinx-argparse
[8]: https://www.python.org
[9]: https://golang.org
[10]: https://www.gonum.org
[11]: https://phylohic.readthedocs.io/en/latest/
[12]: https://www.cnrs.fr/en
[13]: https://tel.archives-ouvertes.fr/tel-02161932
