Experiment Configuration file
=============================

This file contains all the needed configuration to run an experiment. It's
used by SnakeMake as well as the scripts. It's a standard `YAML`_ formatted
file.

.. warning:: All path shall end with a `/`

Here is the content of :file:`config/sample.yaml` that shows all available
configuration keys with a lot of comments. In this sample, we have three
species (suitably named `Species1`, `Species2` and `Species3`), and two Hi-C
resolutions (10kb and 20kb).


.. include:: ../../config/sample.yaml
   :literal:


.. _`YAML`: https://yaml.org
