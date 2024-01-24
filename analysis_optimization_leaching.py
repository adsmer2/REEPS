import numpy as np 
import scipy as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import qsdsan as qs
import biosteam as bst
import os

from systems import create_system
from model_leaching import *

def analysis_optimization_leaching(system, fununit, xdata, ydata, f1data, f2data, indicator, nsamples, figures_path): # set system and functional unit
    # arrays for the full range of each leaching parameter equally spaced
    time = np.linspace(25, 205, nsamples)
    temp = np.linspace(20, 70, nsamples)
    conc = np.linspace(1, 9.5, nsamples)
    solventRatio = np.linspace(2, 7, nsamples)

    # base case/optimal values from Liang et al - 2017
    ftime = 200
    ftemp = 47
    fconc = 3
    fsolventRatio = 2.5

    samples = pd.DataFrame(columns = ['time', 'temp', 'conc', 'solventRatio'])
    
    if ydata == 'time':
        y=time
        ylabel = 'Leaching Time (min)'
        samples['time'] = y
    elif ydata == 'temp':
        y = temp
        ylabel = 'Leaching Temperature (deg C)'
        samples['temp'] = y
    elif ydata == 'conc':
        y = conc
        ylabel = 'Acid Concentration (wt %)'
        samples['conc'] = y
    elif ydata =='solventRatio':
        y = solventRatio
        ylabel = 'Liquid/Solid Ratio'
        samples['solventRatio'] = y
    else:
        raise RuntimeError('"ydata" is not in an acceptable form')

    if f1data == 'time':
        f1=[ftime]*len(y)
        samples['time'] = f1
    elif f1data == 'temp':
        f1 = [ftemp]*len(y)
        samples['temp'] = f1
    elif f1data == 'conc':
        f1 = [fconc]*len(y)
        samples['conc'] = f1
    elif f1data =='solventRatio':
        f1 = [fsolventRatio]*len(y)
        samples['solventRatio'] = f1
    else:
        raise RuntimeError('"f1data" is not in an acceptable form')

    if f2data == 'time':
        f2=[ftime]*len(y)
        samples['time'] = f2
    elif f2data == 'temp':
        f2 = [ftemp]*len(y)
        samples['temp'] = f2
    elif f2data == 'conc':
        f2 = [fconc]*len(y)
        samples['conc'] = f2
    elif f2data =='solventRatio':
        f2 = [fsolventRatio]*len(y)
        samples['solventRatio'] = f2
    else:
        raise RuntimeError('"f2data" is not in an acceptable form')

    if indicator == 'NPV':
        ind_slice = ('TEA', 'NPV15 [MM USD/kg]')
        ind_axis_label = 'NPV (MM USD)'
        ind_cmap = 'inferno'
    elif indicator == 'GWP':
        ind_slice = ('LCA', 'Global Warming [kg CO2-Eq/kg PG]')
        ind_axis_label = f'Global Warming (kg CO2 eq/kg {fununit})'
        ind_cmap = 'inferno_r'
    else:
        raise RuntimeError('"indicator" is not in an acceptable form. Please set as either "NPV" or GWP"')

    model = create_model_leaching(system, fununit) # order of model parameter output [time, temp, conc, solventRatio]
    result = []
    i = 0
    if xdata == 'time':
        x=time
        xlabel = 'Leaching Time (min)'
        while i < len(y):
            xi = [x[i] for n in range(len(y))]
            samples['time'] = xi
            model.load_samples(samples.to_numpy())
            model.evaluate()
            result.extend(model.table.loc[:,ind_slice].values.tolist())
            i += 1
    elif xdata == 'temp':
        x = temp
        xlabel = 'Leaching Temperature (deg C)'
        while i < len(y):
            xi = [x[i] for n in range(len(y))]
            samples['temp'] = xi
            model.load_samples(samples.to_numpy())
            model.evaluate()
            result.extend(model.table.loc[:,ind_slice].values.tolist())
            i += 1
    elif xdata == 'conc':
        x = conc
        xlabel = 'Acid Concentration (wt %)'
        while i < len(y):
            xi = [x[i] for n in range(len(y))]
            samples['conc'] = xi
            model.load_samples(samples.to_numpy())
            model.evaluate()
            result.extend(model.table.loc[:,ind_slice].values.tolist())
            i += 1
    elif xdata =='solventRatio':
        x = solventRatio
        xlabel = 'Liquid/Solid Ratio'
        while i < len(y):
            xi = [x[i] for n in range(len(y))]
            samples['solventRatio'] = xi
            model.load_samples(samples.to_numpy())
            model.evaluate()
            result.extend(model.table.loc[:,ind_slice].values.tolist())
            i += 1
    else:
        raise RuntimeError('"xdata" is not in an acceptable form')

    X, Y = np.meshgrid(x, y)
    Z = np.array(result).reshape(len(y),len(y)).T

    plt.style.use('default')
    aspect_ratio_LtoW = 2880/3840
    cm_to_in = 1/2.54  # centimeters in inches
    width_one_col = 8.3 # cm. Width for a one column figure
    width_two_col = 17.1 # cm. Width for a two column figure
    max_length = 23.3 # cm. The maximum lenght a figure can be
    fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))

    cp = plt.contourf(X,Y,Z, 20, cmap=ind_cmap)
    colorbar = fig.colorbar(cp, label=ind_axis_label)
    # Limit the number of ticks
    num_ticks = 5
    colorbar.locator = matplotlib.ticker.MaxNLocator(nbins=num_ticks)
    # Update the colorbar with the new locator
    colorbar.update_ticks()

    ax.set(xlabel=xlabel, ylabel=ylabel)
    fig.tight_layout()

    fig.savefig(os.path.join(figures_path, f'leaching_{xdata}_{ydata}_{indicator}_{fununit}.tiff'), dpi=600)

    qs.Model._reset_system(model)
