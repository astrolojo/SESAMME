############
File Formats
############


Input Spectra
=============

``SESAMME`` assumes that your integrated-light (IL) star cluster spectrum has already undergone some amount of preprocessing before comparison with the SSP model cube (see the next section `The SSP Model Cube`_). The most important pre-processing item is that your spectrum **must** be in the rest frame. Other recommended preparatory steps include:

* Resampling to the wavelength grid of your SSP models, if necessary. 
* Correction for foreground (Milky Way) extinction.
* Modeling and removal of HI absorption profile (if using observations that cover Ly\ :math:`\alpha`).



The SSP Model Cube
==================

While ``SESAMME`` can accept virtually any set of model SSPs (e.g., those from `BPASS <https://bpass.auckland.ac.nz/index.html>`_ or `Starburst99 <https://www.stsci.edu/science/starburst99/docs/default.htm>`_) to characterize an IL spectrum, it expects those SSPs to be read in a particular format.



