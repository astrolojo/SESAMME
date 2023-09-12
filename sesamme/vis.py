### Manipulating arrays
import numpy as np

### Data visualization
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import corner
from IPython.display import display, Math

plt.rcParams.update({'font.size': 12, 'font.weight':'semibold'})

### Referencing other parts of SESAMME
from sesamme.mcmc import ndim
from sesamme.models import get_model, nebular_continuum, apply_ext_law, use_ext_law


#########################

def print_stats(flat_samples):
    """
    Computes 16th, 50th, and 84th percentile values for the 4 variables and prints them in a friendly way.
    
    Inputs:
    - Flattened MCMC chain flat_samples.
    
    Outputs:
    - Display of formatted Math object.
    """
    global ndim

    labels = ["log(age/yr)", r"log(Z/Z$_{\odot}$)", "E(B-V)", "log(A)"]

    for i in range(ndim):
        percentiles = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(percentiles)
        txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
        txt = txt.format(percentiles[1], q[0], q[1], labels[i])
        display(Math(txt))

#########################

def plot_samples(x, y, windowlist, flat_samples, add_nebular = True, plot_median = False, median_params = [7., -2., 0.0, 0.0], plot_random_draws = True, title = None, savefile = False, savefile_name = "example_fit"):
    """
    A plotting function for examining the goodness-of-fit for models in the final sampler object after the MCMC run.
    
    Inputs:
    - Wavelength array x.
    - Flux array y.
    - List of masked regions windowlist.
    - Flattened MCMC chain flat_samples.
    - Boolean add_nebular, optional. Determines whether to include a nebular continuum component when 
    plotting models. Default is True.
    - Boolean plot_median, optional. Determines whether to plot a single model (usually some sort of best-fit).
    Default is False.
    - List-like or array-like median_params, optional. 
    - Boolean plot_random_draws, optional. Default is True.
    - String title to give a title to the plot, optional. Default is None.
    - Boolean savefile, optional. Default is False.
    - String savefile_name, optional. Specify output file name if savefile set to True. 
    """
    
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10,6), gridspec_kw={'height_ratios': [3, 1]})

    ### Plot the data and masked regions
    ax[0].step(x, y, color='black', lw = 1, label = "Data", zorder=1)

    ### Mark intervals that were masked during fitting
    ### Could probably done a smarter way using the mask object itself instead of windowlist
    for k in range(len(windowlist)):
        maskout = np.where((x >= windowlist[k][0]) & (x <= windowlist[k][1]))
        ax[0].fill_between(x[maskout], 0, 1e6, alpha = 0.5, color='lightgrey', lw=0)
        ax[1].fill_between(x[maskout], -10, 10, alpha = 0.5, color='lightgrey', lw=0)

    ### Optionally plot an individual model (typically a "best fit")
    if plot_median == True:
        ### Uses the same value of add_nebular as given above
        total_model = get_model(median_params, x, modelcube, iontable, add_nebular)
        ax[0].step(x, total_model, color = 'royalblue', lw=1, label = 'Optimal Model', zorder=500)
    
    ### Formatting the upper panel and setting the plot title
    ax[0].set_xlim(np.min(x[mask]-20.), np.max(x[mask]+20.))
    ax[1].set_xlabel(r"Wavelength ($\AA$)")
    
    ylim_low, ylim_high = 0.1 * np.median(y[y != 0.]), 2 * np.median(y[y != 0.])
    ax[0].set_ylim(0, ylim_high)
    ax[0].set_ylabel(r"L$_{\odot}$ $\AA^{-1}$")
    
    if title is not None:
        ax[0].set_title(title, weight='semibold')


    ######

    ### Optionally plot 50 random draws from the final walker ensemble
    if plot_random_draws == True:
        np.random.seed(99)
        for i in np.random.randint(0, len(flat_samples), 50):
            t, z, ebv, amp = flat_samples[i]
            nearest_t, nearest_z = _nearest_age(t)[0], _nearest_metallicity(z)[0]

            draw_y_model = 10**(amp) * modelcube[nearest_z].data[nearest_t]
            
            if add_nebular == True:
                draw_nebcont = nebular_continuum(x, [t, z, ebv, amp], iontable)
            else:
                draw_nebcont = 0
            total_model = draw_y_model + draw_nebcont
            total_model = apply_ext_law(x = x, ebv = ebv, y_model = total_model)

            ax[0].step(x, total_model, alpha=  0.15, lw=1, ls=':', color='teal', zorder=0)
            ax[1].step(x[mask], (y[mask] - total_model[mask])/y[mask], ls=":", alpha=  0.05, lw=1, color='teal', zorder=1)

        
        ### Add label for random draws
        handles, labels = ax[0].get_legend_handles_labels()
        legend_resid = Line2D([0], [0], label='Random Draw from PDF', ls = ":", alpha =0.5, lw=1, color='teal')
        handles.extend([legend_resid])

        ax[0].legend(loc='best', handles=handles)
        
    else:
        ax[0].legend(loc='best')

    
    ### Lower panel formatting
    ax[1].axhline(0, color='black')
    ax[1].set_ylabel("Residuals")
    ax[1].set_ylim(-0.7,0.7)
    
    plt.tight_layout()
    
    
    ### Save the file as a PDF?
    if savefile:
        plt.savefig(savefile_name+".pdf", bbox_inches='tight', overwrite=True)
    plt.show()
