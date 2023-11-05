#####################
Extinction in SESAMME
#####################


Eight extinction curves are currently implemented for ``SESAMME``. All are implemented using either the ``extinction`` package by K. Barbary (`<https://extinction.readthedocs.io/en/latest/index.html>`_) or the ``dust_extinction`` package by K. Gordon (`<https://dust-extinction.readthedocs.io/en/stable/index.html>`_).

+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| Model          | :math:`R_V` | Galaxy Type   |  Reference                                                          | Python source           |
+================+=============+===============+=====================================================================+=========================+
| CCM            | 3.1         | MW            | `<https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract>`_ | ``extinction``          |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| Fitzpatrick99  | 3.1         | MW            | `<https://ui.adsabs.harvard.edu/abs/1999PASP..111...63F/abstract>`_ | ``extinction``          |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| ODonnell       | 3.1         | MW            | `<https://ui.adsabs.harvard.edu/abs/1994ApJ...422..158O/abstract>`_ | ``extinction``          |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| FitzMassa07    | 3.1         | MW            | `<https://ui.adsabs.harvard.edu/abs/2007ApJ...663..320F/abstract>`_ | ``extinction``          |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| Gordon23       | 3.1         | MW            | `<https://ui.adsabs.harvard.edu/abs/2023ApJ...950...86G/abstract>`_ | ``dust_extinction``     |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| Calzetti       | 4.05        | Starburst     | `<https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_ | ``extinction``          |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| SMC            | 2.74        | SMC           | `<https://ui.adsabs.harvard.edu/abs/2003ApJ...594..279G/abstract>`_ | ``dust_extinction``     |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+
| LMC            | 3.41        | LMC           | `<https://ui.adsabs.harvard.edu/abs/2003ApJ...594..279G/abstract>`_ | ``dust_extinction``     |
+----------------+-------------+---------------+---------------------------------------------------------------------+-------------------------+



Parameterizing Extinction
=========================

Although the ``SESAMME`` acronym lists extinction as one of the four variables it models, it is more accurate to say that it samples the line-of-sight *reddening* |EBV|, which can then be translated to an extinction or attenuation |AV| after the fact. This conversion from |EBV| to |AV| depends on the chosen extinction curve and hence the assumed value of |RV|.

For details on how ``SESAMME`` applies reddening values to spectra, see the :ref:`model_gen` section of the API.

.. |EBV| replace:: :math:`E(B-V)`
.. |AV| replace:: :math:`A_V`
.. |RV| replace:: :math:`R_V`


