### File and unit handling tools
from astropy.io import fits
import astropy.units as u
from astropy.table import Table

### Dust models
import extinction, dust_extinction
from dust_extinction.parameter_averages import G23
from dust_extinction.averages import G03_LMCAvg, G03_SMCBar

### Manipulating arrays
import numpy as np
import emcee

### Data visualization
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

import corner
from IPython.display import display, Math

plt.rcParams.update({'font.size': 12, 'font.weight': 'semibold', 'axes.labelsize': 12, 
                     'axes.labelweight':'semibold'})

import warnings



# Set the number and dimensionality of the emcee walker ensemble
# ndim MUST be set to 4
nwalkers, ndim = 128, 4
nsteps = 10000


# Set the initial positions of the walker ensemble
initial_pos = [7., -2.0, 0.2, -2.0] + ([0.1, 0.1, 0.1, 0.1] * np.random.randn(nwalkers, ndim))