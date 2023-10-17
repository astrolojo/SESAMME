### File and unit handling tools
from astropy.table import Table
import warnings

### Manipulating arrays
import numpy as np
import emcee




### Referencing other parts of SESAMME
import sesamme.models as models



# Set the number and dimensionality of the emcee walker ensemble
# ndim MUST be set to 4
nwalkers, ndim = 128, 4
nsteps = 10000


# Set the initial positions of the walker ensemble
initial_pos = [7., -2.0, 0.2, -2.0] + ([0.1, 0.1, 0.1, 0.1] * np.random.randn(nwalkers, ndim))

#########################

def set_chain_size(m):
    """Updates the number of steps (chain length) of the MCMC run. Default value is 10000.
    
    Parameters
    ----------
    m : int
        Number of steps in chain
    """
    global nsteps
    
    nsteps = m

#########################

def set_walker_size(m):
    """Updates the size of the walker ensemble. Default value of nwalkers = 128.
    
    Parameters
    ----------
    m : int
        Number of walkers in emcee ensemble
    """
    global nwalkers
    
    nwalkers = m
    

#########################

def set_initial_positions(centers):
    """Updates the initial positions of the walker ensemble.

    Creates a ball of walker positions with a small spread around specified values. 
    
    Parameters
    ----------
    center : list or np.ndarray
        Must contain, in order, a value for log(age), log(Z), E(B-V), and log(ampl).
    """
    global initial_pos
    
    age_cen, met_cen, ebv_cen, amp_cen = centers[0], centers[1], centers[2], centers[3]
    
    initial_pos = [age_cen, met_cen, ebv_cen, amp_cen,] + ([0.1, 0.1, 0.1, 0.1] * np.random.randn(nwalkers, ndim))


#########################

prior_dict = dict.fromkeys(['age','met', 'ebv', 'amp'])

prior_lowbounds = [6.0, -3.0, 0.01, -20.]
prior_highbounds = [7.5, -1.5, 1.0, 1.0]

def set_prior_bounds(prior_dict, prior_lowbounds, prior_highbounds):
    """Set boundaries on the priors for an MCMC run.
 
    Parameters
    ----------
    prior_dict : dict
        Global dictionary that specifies the boundaries of the four parameters. Initially empty.
    prior_lowbounds : list or np.ndarray
        Gives lower boundaries for the priors. Must contain, in order, a value for log(age), log(Z), E(B-V), and log(ampl).
    prior_highbounds : list or np.ndarray
        Gives upper boundaries for the priors. Must contain, in order, a value for log(age), log(Z), E(B-V), and log(ampl).
    
    Raises
    ------
    ValueError
        If boundaries are not increasing from the first column (expected lower) to the second (higher).
    """

    prior_dict['age'] = [prior_lowbounds[0], prior_highbounds[0]]
    prior_dict['met'] = [prior_lowbounds[1], prior_highbounds[1]]
    prior_dict['ebv'] = [prior_lowbounds[2], prior_highbounds[2]]
    prior_dict['amp'] = [prior_lowbounds[3], prior_highbounds[3]]
    
    _check_prior_bounds(prior_dict)

#########################

def _check_prior_bounds(prior_dict):
    """
    Sanity check that the bounds on the flat priors are laid out sensibly.
    
    Parameters
    ----------
    prior_dict : dict
        Global dictionary that specifies the boundaries of the four parameters.
    
    Raises
    ------
    - ValueError
        If boundaries are not increasing from the first column (expected lower) to the second (higher)
    """
    
    # Ensure that the first value in the prior-boundary pair is lower than the second
    for key, value in prior_dict.items():
        if value[0] > value[1]:
            raise ValueError("Prior boundaries are out of order for variable "+"\'"+key+"\'")
            
    # Check that the prior_dict make physical sense, and throw a warning if they don't
    for key in ['age', 'met', 'ebv']:
        if (key == 'age') & ((prior_dict[key][0] < 5.) | (prior_dict[key][1] < 5.) | (prior_dict[key][0] > 11.) | (prior_dict[key][1] > 11.)):
            warnings.warn("\n\n Your age prior may extend to unphysically young and/or old values. Double check before proceeding with SESAMME.")
        if (key == 'met') & ((prior_dict[key][0] < -5) | (prior_dict[key][1] < -5) | (prior_dict[key][0] > -1.3) | (prior_dict[key][1] > -1.3)):
            warnings.warn("\n\n Your metallicity prior may extend to unphysical values. Double check before proceeding with SESAMME.")
        if (key == 'ebv') & ((prior_dict[key][0] < 0) | (prior_dict[key][1] < 0) | (prior_dict[key][0] > 100) | (prior_dict[key][1] > 100)):
            warnings.warn("\n\n Your E(B-V) prior may extend to unphysical values. Double check before proceeding with SESAMME.")


#########################

def log_prior(theta):
    """
    Calculate a flat prior probability within the bounds prescribed below.
    
    Parameters
    ----------
    theta : list or np.ndarray
        Array containing, in order, a log(age) logt, a log(metallicity) logZ, an E(B-V) value ebv, and an amplitude ampl.
    

    Returns
    -------
    float
        0 if parameters in ``theta`` are within allowed ranges, or -inf otherwise
    """
    
    # Ensure the priors make sense...
    _check_prior_bounds(prior_dict)
    
    # ...Then use them to set the flat prior probability
    
    logt, logZ, ebv, ampl = theta
    
    if prior_dict['age'][0] <= logt <= prior_dict['age'][1] and \
    prior_dict['met'][0] <= logZ <= prior_dict['met'][1] and \
    prior_dict['ebv'][0] < ebv < prior_dict['ebv'][1] and \
    prior_dict['amp'][0] < ampl < prior_dict['amp'][1] :
        return 0.0
    else:
        return -np.inf

#########################

def log_likelihood(y, yerr, y_model, mask):
    """
    Evaluate the likelihood function by comparing the data and chosen model at every wavelength bin. 
    
    Parameters
    ----------
    y : array-like
        Flux array
    yerr : array-like
        Flux uncertainty array
    y_model : array-like
        Model flux array
    mask : array-like
        Marks which wavelength bins to ignore during fitting
    theta : list or np.ndarray
        Array containing, in order, a log(age) logt, a log(metallicity) logZ, an E(B-V) value ebv, and an amplitude ampl.
    

    Returns
    -------
    resid : float
        -1 * log(likelihood), where log(n) means the natural logarithm
    """  
    
    masked_spec = np.array( (y[mask] - y_model[mask]) / yerr[mask] )
    masked_err = np.array( np.sqrt(2) * np.sqrt(np.pi) * yerr[mask])
    
    resid = -0.5 * (np.dot(masked_spec, masked_spec) + np.log(np.dot(masked_err, masked_err)))

    return resid

#########################

def log_posterior(theta, x, y, yerr, model_cube, ion_table, mask, add_nebular):
    """
    Calculates the (log of the) posterior probability as log(Ppos) = log(Pprior) + log(likelihood).
 
    Parameters
    ----------
    theta : list or np.ndarray
        Array containing, in order, a log(age) logt, a log(metallicity) logZ, an E(B-V) value ebv, and an amplitude log(ampl).
    x : array-like
        Wavelength array
    y : array-like
        Flux array
    yerr : array-like
        Flux uncertainty array
    model_cube : FITS
        Multi-extension FITS cube containing SSP models
    ion_table : astropy Table
        Table object containing ionizing fluxes per SSP
    mask : array-like
        Marks which wavelength bins to ignore during fitting
    add_nebular : Boolean
        Determines whether to add nebular continuum emission to a model

    Returns
    -------
    lp_ll : float
        Sum of ln(prior) and ln(likelihood)
    """
    y_model = models.get_model(theta, x, model_cube, ion_table, add_nebular)
    
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    
    ll = log_likelihood(y, yerr, y_model, mask)
    
    return lp + ll

#########################

def run_sesamme(filename, runname, x, y, yerr, model_cube, ion_table, mask, add_nebular=True):
    """
    Initiate an MCMC procedure and write the results to a file/extension name.
    

    Parameters
    ----------
    filename : str
        Name of \*.H5 file to write results to
    runname : str
        Extension/run name within the file; can store multiple runs in a single \*.H5 file
    x : array-like
        Wavelength array
    y : array-like
        Flux array
    yerr : array-like
        Flux uncertainty array
    model_cube : FITS
        Multi-extension FITS cube containing SSP models
    ion_table : astropy Table
        Table object containing ionizing fluxes per SSP
    mask : array-like
        Marks which wavelength bins to ignore during fitting
    add_nebular : Boolean
        Determines whether to add nebular continuum emission to a model
    """

    global nwalkers, nsteps, ndim
    global initial_pos

    # Confirm with the user that they want to proceed with the current settings
    print("Active extinction law = "+models.use_ext_law + "; Ensemble size = "+str(nwalkers) + 
         "; Chain length = "+str(nsteps))
    yesno = input("Begin a SESAMME run with these parameters? (y/n)  ")
    
    if yesno in ['n', 'no', 'N', 'NO']:
        print("Run canceled.")
        return
    
    elif yesno in ['y', 'Y', 'yes', 'YES']:
        # Set up new HDF backend and run the EnsembleSampler in the usual emcee fashion
        backend = emcee.backends.HDFBackend(filename, name = runname)
        backend.reset(nwalkers, ndim)
        
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_posterior, backend = backend, args=(x, y, yerr, 
                                                  model_cube, ion_table, mask, add_nebular)
        )
        
        sampler.run_mcmc(initial_pos, nsteps, progress=True);
        
        print(
            "Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction))
        )
        
    else:
        print("\n Unknown response. Please use y / yes / Y / YES to begin or n / no / N / NO to change your mind.")
