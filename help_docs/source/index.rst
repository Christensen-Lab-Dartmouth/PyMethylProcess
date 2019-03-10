.. PyMethylProcess documentation master file, created by
   sphinx-quickstart on Sun Dec  9 18:45:51 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyMethylProcess's documentation!
=====================================

.. image:: pymethylprocess_overview.jpeg
   :width: 800px
   :height: 600px
   :alt: PyMethylProcess
   :align: center

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. image:: pipeline-download.jpeg
   :width: 800px
   :height: 600px
   :scale: 60%
   :alt: Download
   :align: center

.. image:: pipeline-format.jpeg
   :width: 800px
   :height: 600px
   :scale: 60%
   :alt: Format
   :align: center

.. image:: pipeline-preprocess.jpeg
   :width: 800px
   :height: 600px
   :scale: 60%
   :alt: Preprocess
   :align: center

.. image:: pipeline-visualize.jpeg
   :width: 800px
   :height: 600px
   :scale: 60%
   :alt: Visualize
   :align: center

.. image:: pipeline-train-test-split.jpeg
   :width: 800px
   :height: 600px
   :scale: 60%
   :alt: TrainTestSplit
   :align: center

.. argparse::
   :module: basic_installer
   :func: main
   :prog: pymethyl-basic-install

.. click:: installer:install
   :prog: pymethyl-install
   :show-nested:

.. click:: visualizations:visualize
   :prog: pymethyl-visualize
   :show-nested:

.. click:: preprocess:preprocess
   :prog: pymethyl-preprocess
   :show-nested:

.. click:: utils:util
   :prog: pymethyl-utils
   :show-nested:

.. argparse::
   :module: run_random_forest
   :func: main
   :prog: pymethyl-basic-ml

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
