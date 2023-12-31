### File and unit handling tools
from astropy.io import fits
from astropy.table import Table
import astropy.units as u

### Manipulating arrays
import numpy as np
from scipy import interpolate

### Dust models
import extinction, dust_extinction
from dust_extinction.parameter_averages import G23
from dust_extinction.averages import G03_LMCAvg, G03_SMCBar

use_ext_law = 'CCM'


def load_ssp_cube(file_name):
    """
    Loads in the SSP model cube.
    Implicitly alters the model grid to be sampled with emcee through _ingest_model_grid().

    Parameters
    ----------
    file_name : str
        File path and name

    Returns
    -------
    model_cube : FITS
        Multi-extension FITS cube containing SSP models
    """

    model_cube = fits.open(file_name)
    
    _ingest_model_grid(model_cube)
    
    return model_cube

def load_ionization_table(file_name):
    """
    Loads in the table of ionizing photon outputs per SSP associated with a model cube.

    Parameters
    ----------
    file_name : str
        File path and name

    Returns
    -------
    ion_table : astropy Table
        Table object containing ionizing fluxes per SSP
    """

    ion_table = Table.read(file_name, format='ascii', header_start=0)
    
    return ion_table

#########################


def _interpret_metallicity_keys(model_cube):
    """
    Interprets the FITS extension names of the model cube to create the numerical metallicity grid. 

    Parameters
    ----------
    model_cube : FITS
        Multi-extension FITS cube containing SSP models

    Returns
    -------
    metal_values : list
        List of metallicity values describing the model cube
    """

    metal_keys = [model_cube[k].header['EXTNAME'] for k in range(1, len(model_cube))]
    metal_values = []
    
    for key in metal_keys:
    
        if 'em' in key:
            _, power = key.split('em')
            met_val = float("1e-" + power)


        else:
            _, value = key.split('Z')
            met_val = float("." + value)
            
        metal_values.append(met_val)

    return metal_values

#########################


def _ingest_model_grid(model_cube):
    """
    Creates the metallicity and age dictionaries that are needed to translate between MCMC samples and discrete age+Z combos that exist in the model cube.


    Parameters
    ----------
    model_cube : FITS
        Multi-extension FITS cube containing SSP models
    """

    global metal_dict, age_dict
    
    metal_dict = {}
    age_dict = {}
    
    example_table = model_cube[1].data
    
    for a in example_table.names[1:]:
        age_dict[a] = float(a)
        
    metal_headers = [model_cube[k].header['EXTNAME'] for k in range(1, len(model_cube))]
    metal_values = _interpret_metallicity_keys(model_cube)
    metal_dict = dict(zip(metal_headers, metal_values))

#########################

### Round to the nearest log(age) and log(metallicity) that is precomputed in the input SSP suite.

def _nearest_age(logt):
    """
    Rounds the input value of log(age) to the nearest value in the current model grid.


    Parameters
    ----------
    logt : float
        Log of cluster age in yr

    Returns
    -------
    best_t_key : str
        Dict key for nearest age
    best_t : float
        Dict value for nearest age
    """
    global age_dict

    best_t_key, best_t = min(age_dict.items(), key=lambda x: abs(logt - x[1]))
    return best_t_key, best_t

def _nearest_metallicity(logZ):
    """
    Rounds the input value of log(Z) to the nearest value in the current model grid.


    Parameters
    ----------
    logZ : float
        Log of stellar metallicity

    Returns
    -------
    best_Z_key : str
        Dict key for nearest metallicity
    best_Z : float
        Dict value for nearest metallicity
    """

    global metal_dict

    Z = 10**logZ
    best_Z_key, best_Z = min(metal_dict.items(), key=lambda x: abs(Z - x[1]))
    return best_Z_key, best_Z

#########################

def get_mask(windowlist, x):
    """
    Define one or more subsets of wavelength space to ignore when performing the data-model comparison (e.g., sky or geocoronal lines, ISM features).

    Each row in the input ``windowlist`` specifies a wavelength window where weight = 0 during the fitting, and the first and second columns respectively give the approximate lower and upper bounds of the window(s). The actual bounds will be chosen from the nearest values in the wavelength array.

    Parameters
    ----------
    windowlist : np.ndarray
       N x 2 array. 
    x : array-like
        Wavelength array

    Returns
    -------
    mask_array : np.ndarray
        Array of bools (1 = use in fit, 0 = do not use)
    """

    mask_array = np.ones(len(x))
    
    windowmin, windowmax = np.zeros(len(windowlist)), np.zeros(len(windowlist))
    
    for k in range(len(windowlist)):
        idxmin = (np.abs(x - windowlist[k][0])).argmin()
        idxmax = (np.abs(x - windowlist[k][1])).argmin()
        windowmin[k] += x[idxmin]
        windowmax[k] += x[idxmax]
        
        mask_array[np.where( (x >= windowmin[k]) & (x <= windowmax[k]) )] = 0
    
    return np.array(mask_array, dtype='bool')

#########################

def set_ext_law(ext_curve):
    """
    Sets the curve used to extinguish model spectra. 

    Currently implemented options for extinction curves include 5 Milky Way-like curves (CCM, Fitzpatrick99, ODonnell, FitzMassa07, Gordon23), a starburst galaxy curve (Calzetti), an LMC-like curve (LMC), and an SMC-like curve (SMC).

    Parameters
    ----------
    ext_curve : str
        Name of extinction curve. Must match an implemented option in SESAMME's library.

    Raises
    ------
    ValueError
        String ext_curve is not in the list of allowable values
    """ 
    
    global use_ext_law
    
    if ext_curve not in ['CCM', 'Fitzpatrick99', 'ODonnell', 'FitzMassa07', 'Gordon23',
                           'Calzetti', 'SMC', 'LMC']:
        raise ValueError("'"+ext_curve+"'" +" is not a valid choice of extinction law.\n \
Accepted values are 'CCM', 'Fitzpatrick99', 'ODonnell', 'FitzMassa07', 'Gordon23', Calzetti', 'LMC', and 'SMC'.")
    
    use_ext_law = ext_curve
    
    print("Model spectra will now be reddened assuming the", ext_curve, "extinction curve")

#########################

def apply_ext_law(x, ebv, y_model):
    """
    Redden a model spectrum for comparison with the data. 
    
    Parameters
    ----------
    x : array-like
        Wavelength array
    ebv : float
        Value of E(B-V) by which to redden the model
    y_model : np.ndarray
        Flux array of the model spectrum; may be purely stellar or stellar + nebular continuum

    Returns
    -------
    red_model : array-like
        Extinguished SSP or SSP+nebular model

    Raises
    ------
    ValueError
        String use_ext_law is not in the list of allowable values

    """    
    
    if use_ext_law not in ['CCM', 'Fitzpatrick99', 'ODonnell', 'FitzMassa07', 'Gordon23',
                           'Calzetti', 'SMC', 'LMC']:
        raise ValueError("'"+use_ext_law+"'" +" is not a valid choice of extinction law.\n \
Accepted values are 'CCM', 'Fitzpatrick99', 'ODonnell', 'FitzMassa07', 'Gordon23', Calzetti', 'LMC', and 'SMC'.")
        
    elif use_ext_law == "CCM":
        red_model = extinction.apply(extinction.ccm89(x, ebv*3.1, 3.1), y_model)
        
    elif use_ext_law == "Fitzpatrick99":
        red_model = extinction.apply(extinction.fitzpatrick99(x, ebv*3.1, 3.1), y_model)
        
    elif use_ext_law == "ODonnell":
        red_model = extinction.apply(extinction.odonnell94(x, ebv*3.1, 3.1), y_model)
        
    elif use_ext_law == "FitzMassa07":
        red_model = extinction.apply(extinction.fm07(x, ebv*3.1), y_model)
        
    elif use_ext_law == "Calzetti":
        red_model = extinction.apply(extinction.calzetti00(x, ebv*4.05, 4.05), y_model)
    
    elif use_ext_law == 'Gordon23':
        ext_curve = G23(Rv = 3.1)
        red_model = y_model * ext_curve.extinguish(x * u.AA, ebv * ext_curve.Rv)
        
    elif use_ext_law == 'LMC':
        ext_curve = G03_LMCAvg()
        red_model = y_model * ext_curve.extinguish(x * u.AA, ebv * ext_curve.Rv)
    
    elif use_ext_law == 'SMC':
        ext_curve = G03_SMCBar()
        red_model = y_model * ext_curve.extinguish(x * u.AA, ebv * ext_curve.Rv)
        
    return red_model

#########################

def get_model(theta, x, model_cube, ion_table, add_nebular = True):
    """
    Construct a model star cluster spectrum with values of metallicity, age, extinction, and normalization randomly chosen (within bounds).

    Choose the model nearest to the age and metallicity parameters Z and t, rescale by ampl, then redden by the extinction parameter ebv. If add_nebular is set to True, this is added before rescaling and extinction.

    Parameters
    ----------
    theta : list or np.ndarray
        Array containing, in order, a log(age) logt, a log(metallicity) logZ, an E(B-V) value ebv, and an amplitude log(ampl).
    model_cube : FITS
        Multi-extension FITS cube containing SSP models
    ion_table : astropy Table
        Table object containing ionizing fluxes per SSP
    add_nebular : Boolean
        Determines whether to add nebular continuum emission to a model

    Returns
    -------
    red_nearest_model : array-like
        Extinguished and rescaled SSP or SSP+nebular model spectrum with age and metallicity values nearest to the inputs. 
    """
    
    logt, logZ, ebv, ampl = theta
    
    met_key, met = _nearest_metallicity(logZ)
    age_key, age = _nearest_age(logt)
    
    nearest_stellar_model = (10**ampl) * model_cube[met_key].data[age_key]
    specwl = model_cube[met_key].data['WL'].byteswap().newbyteorder()
    
    if add_nebular == True:
        nebcont = nebular_continuum(x, theta, ion_table)
        nearest_model = nearest_stellar_model + nebcont
    else:
        nearest_model = nearest_stellar_model
    
    red_nearest_model = apply_ext_law(x, ebv, nearest_model)        
    
    return red_nearest_model

#########################

gamma = np.array([0.,2.11e-4,5.647,9.35,9.847,10.582,16.101,24.681,26.736,
              24.883,29.979,6.519,8.773,11.545,13.585,6.333,10.444,7.023,
              9.361,7.59,9.35,8.32,9.53,8.87])*1e-40  # Units are 10**-40 ERG CM3 SEC**-1 HZ**-1
nebx = np.array([912.,913.,1300.,1500.,1800.,2200.,2855.,3331.,3421.,3422.,
             3642.,3648.,5700.,7000.,8207.,8209.,14583.,14585.,22787.,22789.,
             32813.,32815.,44680.,44682.])
alpha_B = 2.6e-13    
Qbase = 52

sparse_nebcont = (2.998e18 * gamma * 10**(Qbase)) / (alpha_B * nebx**2)


def nebular_continuum(x, theta, ion_table):
    """
    Python equivalent of the Starburst99 function CONTINUUM, for computing the approximate nebular continuum.
    
    Variable 'gamma' contains emission coeffs for HI (including free-free, bound-free, and 2-photon emission) and HeI (assuming He/H = 0.1). Values from Aller (1984) and Ferland (1980). All computations assume typical Case B conditions (T = 1e4 K,  f_esc = 0).

    Parameters
    ----------
    x : array-like
        Wavelength array
    theta : list or np.ndarray
        Array containing, in order, a log(age) logt, a log(metallicity) logZ, an E(B-V) value ebv, and an amplitude log(ampl).
    ion_table : astropy Table
        Table object containing ionizing fluxes per SSP

    Returns
    -------
    interp_nebcont : np.ndarray
        Approximate rescaled nebular continuum spectrum in Solar luminosities per A.
    """

    logt, logZ, _, ampl = theta
    met_key, met = _nearest_metallicity(logZ)
    age_key, age = _nearest_age(logt)
    
    # Interpolate the sparsely-sampled continuum to the wavelength grid of the model
    interp_function = interpolate.interp1d(nebx, sparse_nebcont, fill_value='extrapolate', )
    interp_nebcont = interp_function(x) / 3.83e33
    
    # Rescale the continuum by the emissivity of the nearest model and by the specified amplitude parameter
    Q_new = ion_table[ion_table['Z'] == met_key][age_key] 
    interp_nebcont = interp_nebcont * (10**(Q_new - Qbase + ampl))  
    
    return interp_nebcont