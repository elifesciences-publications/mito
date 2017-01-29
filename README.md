mito
====

Tools for analysing mitochondrial distribution in fluoresence microscopy images

This code is a little gnarly. In an effort to reduce memory usage it makes
heavy use of iterators (actually generator expressions), which means all
computation is based on a single pass through the iterator. This means the code
is quite nested.


* Free software: MIT license


Installation
------------

```
pip install -U numpy
pip install --process-dependency-links -e git+https://gitlab.erc.monash.edu.au/skeith/mito.git#egg=mito
```

To use the notebook:

```
pip install -U jupyter
```

Features
--------

Currently, the best way to run this analysis is via a `jypter notebook`. The analysis runs specifically on multi-series Leica LIF files and primarily implements the following routines:

* Calculate R90% (images need DAPI, MitoTracker and RSV1-GFP channels).
* Analyse distribution of mitochondria relative to the MTOC (DAPI, MitoTracker,
  RSV1-GFP and gamma-tubulin channels).

Basic Usage
-----------


Imports:

```python

  # For mito analysis functions and routines
  from mito import mito

  # For plotting functions
  from mito import plot
```

Calculate R90%:

```python

  import javabridge
  import bioformats
  import pyome

  from mito import mito

  img_path = "/path/to/img.lif"
  output_dir = "/path/to/output_dir"

  r90 = process_r90(img_path, treatment="Mock", time="8h", output_dir=output_dir)
```

Full examples showing how to run the analysis are in the `mito_analysis.ipynb` notebook in the notebook directory.

```
cd notebooks
jupyter notebooks
```

Note: To run the notebook, you will need `jupyter notebook` as described in
[Installation][#Installation]. Images for this study are available on request. You will need to update the paths to the location of the images provided.

Credits
---------

This package was created with [Cookiecutter][cc] and the [audreyr/cookiecutter-pypackage][aud] project template.

[cc]: https://github.com/audreyr/cookiecutter
[aud]: https://github.com/audreyr/cookiecutter-pypackage

