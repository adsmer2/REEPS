import numpy as np 
import scipy as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import qsdsan as qs
import biosteam as bst
import os

from systems import create_system
from model import *

def analysis_indicator_trend(model, p_name, figures_path):
    '''
    Create a line plot that shows the change in an indicator over the uncertainty range of a parameter

    model = any model that has been run previously

    p_name = the name as designated in the model

    figures_path = ospath to print the figure to (as .tiff at 600 dpi)
    '''
    # Indicator Trend Analysis
    # --------------------------
    
    # Define parameter of interest
    ind_parameter = [i.name for i in model._parameters].index(p_name)
    p_units = model._parameters[ind_parameter].units
    
    # Define range of parameter
    lower = model._parameters[ind_parameter].distribution.lower[0]
    upper = model._parameters[ind_parameter].distribution.upper[0]
    num_bins = 65
    
    # Define the bin width and make a list of bin centers
    bin_width = (upper-lower)/num_bins
    bin_centers = np.linspace(lower + bin_width/2, upper - bin_width/2, num=int((upper - lower) / bin_width))
    mean_output_parameter_binned = []
    std_output_parameter_binned = []

    # Gather the mean and standard deviation for each bin from the model results
    for bin_center in bin_centers:
        bin_start = bin_center - bin_width/2
        bin_end = bin_center + bin_width/2
        bin_indices = np.logical_and(model.table.loc[:, (model._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values >= bin_start, 
                                    model.table.loc[:, (model._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values < bin_end)
        bin_data = model.table.loc[:,('TEA', 'IRR [%]')].values[bin_indices]
        mean_output_parameter_binned.append(np.mean(bin_data))
        std_output_parameter_binned.append(np.std(bin_data))

    # Smooth the data using cubic spline interpolation
    smooth_x_values = np.linspace(lower, upper, num=1000)  # Increase the number of points for smooth curves
    mean_interpolator = np.interp(smooth_x_values, bin_centers, mean_output_parameter_binned)
    std_interpolator = np.interp(smooth_x_values, bin_centers, std_output_parameter_binned)

    # Create the line plot with uncertainty using plot and fill_between
    fig, ax = plt.subplots()
    ax.plot(smooth_x_values, mean_interpolator, color='#60c1cf', label='Mean IRR')
    ax.fill_between(smooth_x_values, mean_interpolator - 2*std_interpolator, mean_interpolator + 2*std_interpolator, color='#60c1cf', alpha=0.3, label='Confidence Interval (95%)')
    # Add labels and title
    ax.set(xlabel=f'{p_name} ({p_units})', ylabel='IRR (fraction)', xlim=(lower,upper), ylim=(0,0.4))
    # Show the legend
    ax.legend(loc='upper left')
    # Save the plot
    fig.savefig(os.path.join(figures_path, f'IRR_Trend_{p_name}.tiff'), dpi=600)