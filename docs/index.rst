#######
SESAMME
#######

``SESAMME`` is a Python package for Simultaneous Estimates of Star-cluster Age, Metallicity, Mass, and Extinction. It uses Bayesian methods to model the integrated-light spectrum of a star cluster, assuming a user-specified set of simple stellar population (SSP) models. 

This package uses the ``emcee`` framework for its Markov Chain Monte Carlo computations. For additional details on the ``SESAMME`` algorithm and a look at its first science applications, see `Jones, L. H. et al. 2023`.



User Documentation
==================

.. note::
	This documentation is under active development

.. toctree::
   :maxdepth: 2

   File Formats <file_formats.rst>
   Extinction in SESAMME <extinction_curves.rst>

Installation
============

.. toctree::
   :maxdepth: 2

   How to install <install.rst>

Source code can be retrieved from `@astrolojo on Github <https://github.com/astrolojo/SESAMME>`_


Reporting Issues
================

If you have found a bug in ``SESAMME`` please report it by creating a
new issue on the ``SESAMME`` `GitHub issue tracker
<https://github.com/astrolojo/SESAMME/issues>`_.

Please include an example that demonstrates the issue sufficiently so that the
developers can reproduce and fix the problem.


.. _reference_API:

Reference API
=============

.. toctree::
   :maxdepth: 2

   API <api.rst>
