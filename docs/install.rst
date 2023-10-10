##############
How to Install
##############

From source
===========

``SESAMME`` can be installed from the source code in the normal
Python fashion after downloading it from the git repo::

    pip install -e .

Using pip
=========

``SESAMME`` can also be installed using pip::

    # from the main branch of the repository
    pip install git+https://github.com/astrolojo/SESAMME.git

If you get an error that is includes `SSLError(SSLCertVerificationError`, the
following command may work to remove this error::

    pip install --trusted-host=pypi.org --trusted-host=files.pythonhosted.org git+https://github.com/astrolojo/SESAMME.git
