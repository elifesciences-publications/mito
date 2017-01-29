===============================
mito
===============================


.. https://img.shields.io/pypi/v/mito.svg
        :target: https://pypi.python.org/pypi/mito

.. https://img.shields.io/travis/keithschulze/mito.svg
        :target: https://travis-ci.org/keithschulze/mito

.. https://readthedocs.org/projects/mito/badge/?version=latest
        :target: https://mito.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Tools for analysing mitochondrial distribution in fluoresence microscopy images

This code is a little gnarly. In an effort to reduce memory usage it makes
heavy use of iterators (actually generator expressions), which means all
computation is based on a single pass through the iterator. This means the code
is quite nested.


* Free software: MIT license


Installation
------------

..code-block:: bash

  pip install -U numpy
  pip install --process-dependency-links -e https://gitlab.erc.monash.edu.au/skeith/mito.git


To use the notebook:

..code-block:: bash

  pip install -U jupyter


Features
--------

Currently, the best way to run this analysis is via a `jypter notebook`. The analysis runs specifically on multi-series Leica LIF files and primarily implements the following routines:

* Calculate R90% (images need DAPI, MitoTracker and RSV1-GFP channels).
* Analyse distribution of mitochondria relative to the MTOC (DAPI, MitoTracker,
  RSV1-GFP and gamma-tubulin channels).

Basic Usage
-----------

Full examples showing how to run the analysis are in the `mito_analysis.ipynb` notebook in the notebook directory.

Imports:

..code-block:: python

  # For mito analysis functions and routines
  from mito import mito

  # For plotting functions
  from mito import plot


Calculate R90%:

..code-block:: python

  import javabridge
  import bioformats
  import pyome

  from mito import mito

  img_path = "/path/to/img.lif"
  output_dir = "/path/to/output_dir"

  r90 = process_r90(img_path, treatment="Mock", time="8h", output_dir=output_dir)


Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

