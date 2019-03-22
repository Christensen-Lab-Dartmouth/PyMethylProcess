Welcome to PyMethylProcess's documentation!
===========================================

https://github.com/Christensen-Lab-Dartmouth/PyMethylProcess

To get started, download pymethylprocess using Docker (joshualevy44/pymethylprocess) or PIP (pymethylprocess) and run pymethyl-install_r_dependencies.

There is both an API and CLI available for use. Examples for CLI usage can be found in ./example_scripts.

# todo: add images, add other two CLI, fix general machine learning, document all classes, publish and reference published help docs in readme

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. image:: .images/pipeline-download.jpeg
  :width: 800px
  :height: 600px
  :scale: 60%
  :alt: Download
  :align: center

.. image:: .images/pipeline-format.jpeg
  :width: 800px
  :height: 600px
  :scale: 60%
  :alt: Format
  :align: center

.. image:: .images/pipeline-preprocess.jpeg
  :width: 800px
  :height: 600px
  :scale: 60%
  :alt: Preprocess
  :align: center

.. image:: .images/pipeline-visualize.jpeg
  :width: 800px
  :height: 600px
  :scale: 60%
  :alt: Visualize
  :align: center

.. image:: .images/pipeline-train-test-split.jpeg
  :width: 800px
  :height: 600px
  :scale: 60%
  :alt: TrainTestSplit
  :align: center

.. automodule:: pymethylprocess.PreProcessDataTypes
  :members:

.. automodule:: pymethylprocess.MethylationDataTypes
  :members:

.. automodule:: pymethylprocess.meffil_functions
  :members:

.. click:: pymethylprocess.installer:install
  :prog: pymethyl-install
  :show-nested:

.. click:: pymethylprocess.visualizations:visualize
  :prog: pymethyl-visualize
  :show-nested:

.. click:: pymethylprocess.preprocess:preprocess
  :prog: pymethyl-preprocess
  :show-nested:

.. click:: pymethylprocess.utils:util
  :prog: pymethyl-utils
  :show-nested:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
