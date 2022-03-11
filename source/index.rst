.. MATTE documentation master file, created by
   sphinx-quickstart on Wed Mar  9 14:06:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MATTE's documentation!
=================================

Welcome to raise issues in github https://github.com/gx-Cai/MATTE

Description
-----------

MATTE (Module Alignment of TranscripTomE) is a python package aiming to
analysis transcriptome from samples with different phenotypes in a
module view. Differiential expression (DE) is commonly used in analysing
transcriptome data. But genes are not work alone, they collaborate.
Network and module based differential methods are developed in recent
years to obtain more information. New problems appears that how to make
sure module or network structure is preserved in all of the phenotypes.
To that end, we proposed MATTE to find the conserved module and diverged
module by treating genes from different phenotypes as individual ones.
By doing so, meaningful markers and modules can be found to better
understand what's really difference between phenotypes.

**Advantages**

In the first place, MATTE merges the data from phenotypes, seeing genes
from different phenotypes as new analyzing unite. By doing so, benefits
got as follows:

1. MATTE considering the information in phenotypes in the preprocessing
   stage, hoping to find more interesting conclusion.
2. MATTE is actually making transcriptome analysis includes the
   relationship between phenotypes, which is of significance in cancer
   or other complex phenotypes.
3. MATTE can deal with more noise thanks to calculation of relative
   different expression (RDE) and ignore some of batch effect.
4. In a module view, “Markers” can be easily transfer to other case but
   not over fits compare to in a gene view.
5. The result of MATTE can be easily analysed.

Install
-------

Install from pip is recommended.

::

   pip install MATTE


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Contents:
.. toctree:: 
   MATTE
   Usergide
