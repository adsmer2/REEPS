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

# =============================================================================
# Target Analysis
# =============================================================================
def analysis_target(fununit, feedPG, REEcontent, num_ind_REEs, num_samples, figures_path):
    # Target 1
    # -------------------
    sys, lca, tea= create_system(fununit=fununit, feedPG=feedPG, REEcontent=REEcontent, num_ind_REEs=num_ind_REEs)
    target_n = 'Adsorbent Capacity (S1)'
    
    model2 = create_model(sys, fununit, 'all', target=target_n) # set system, functional unit, and type of parameter for analysis

    np.random.seed(3221) # setting the seed ensures you get the same sample

    samples = model2.sample(N=num_samples, rule='L')
    model2.load_samples(samples)
    model2.evaluate()

    # Define parameter of interest
    p_name = target_n
    ind_parameter = [i.name for i in model2._parameters].index(p_name)
    p_units = model2._parameters[ind_parameter].units
    
    # Define range of parameter
    lower = model2._parameters[ind_parameter].distribution.lower[0]
    upper = model2._parameters[ind_parameter].distribution.upper[0]
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
        bin_indices = np.logical_and(model2.table.loc[:, (model2._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values >= bin_start, 
                                    model2.table.loc[:, (model2._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values < bin_end) # Gather indices of the parameter values within the bin
        bin_data = model2.table.loc[:,('TEA', 'NPV15 [MM USD/kg]')].values[bin_indices] # using the relevant bin indices, create a list of indicator values.  .loc[:,('TEA', 'IRR [%]')]
        mean_output_parameter_binned.append(np.mean(bin_data)) # take the mean of the indicator values for this bin
        std_output_parameter_binned.append(np.std(bin_data)) # take the SD of the indicator values for this bin
    # Smooth the data using cubic spline interpolation 
    # Note: interp fills the gap between 'lower'/smallest 'bin_center' and 'upper'/higher 'bin_center' leading to a small 'flatline' at either end of plot. Can be minimized with more bins.
    smooth_x_values = np.linspace(lower, upper, num=1000)  # Increase the number of points for smooth curves
    mean_interpolator = np.interp(smooth_x_values, bin_centers, mean_output_parameter_binned)
    std_interpolator = np.interp(smooth_x_values, bin_centers, std_output_parameter_binned)

    # Create the line plot with uncertainty using plot and fill_between
    plt.style.use('default')
    fig2, ax2 = plt.subplots()
    ax2.plot(smooth_x_values, mean_interpolator, color='#f98f60', label='Mean NPV15')
    ax2.fill_between(smooth_x_values, mean_interpolator - 2*std_interpolator, mean_interpolator + 2*std_interpolator, color='#f98f60', alpha=0.3, label='Confidence Interval (95%)')
    # Add labels and title
    ax2.set(xlabel=f'{p_name} ({p_units})', ylabel='NPV15 (MM USD)', xlim=(lower,upper))
    ax2.grid(visible=True)
    ax2.ticklabel_format(axis='x', style='scientific', scilimits=(-3,-2))
    # Show the legend
    ax2.legend(loc='upper left')
    ax2.xaxis.major.formatter._useMathText = True
    fig2.tight_layout()
    # Save the plot
    fig2.savefig(os.path.join(figures_path, f'Target_{p_name}.tiff'), dpi=600)


    # Target 2
    # --------------------------
    sys, lca, tea= create_system(fununit=fununit, feedPG=feedPG, REEcontent=REEcontent, num_ind_REEs=num_ind_REEs)
    target = 'REE Recovery (S1)'
    model = create_model(sys, fununit, 'all', target=target) # set system, functional unit, and type of parameter for analysis

    np.random.seed(3221) # setting the seed ensures you get the same sample

    samples = model.sample(N=num_samples, rule='L')
    model.load_samples(samples)
    model.evaluate()

    # Define parameter of interest
    p_name = target
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
                                    model.table.loc[:, (model._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values < bin_end) # Gather indices of the parameter values within the bin
        bin_data = model.table.loc[:,('TEA', 'NPV15 [MM USD/kg]')].values[bin_indices] # using the relevant bin indices, create a list of indicator values.  .loc[:,('TEA', 'IRR [%]')]
        mean_output_parameter_binned.append(np.mean(bin_data)) # take the mean of the indicator values for this bin
        std_output_parameter_binned.append(np.std(bin_data)) # take the SD of the indicator values for this bin
    # Smooth the data using cubic spline interpolation 
    # Note: interp fills the gap between 'lower'/smallest 'bin_center' and 'upper'/higher 'bin_center' leading to a small 'flatline' at either end of plot. Can be minimized with more bins.
    smooth_x_values = np.linspace(lower, upper, num=1000)  # Increase the number of points for smooth curves
    mean_interpolator = np.interp(smooth_x_values, bin_centers, mean_output_parameter_binned)
    std_interpolator = np.interp(smooth_x_values, bin_centers, std_output_parameter_binned)

    # Create the line plot with uncertainty using plot and fill_between
    plt.style.use('default')
    fig, ax = plt.subplots()
    ax.plot(smooth_x_values, mean_interpolator, color='#f98f60', label='Mean NPV15')
    ax.fill_between(smooth_x_values, mean_interpolator - 2*std_interpolator, mean_interpolator + 2*std_interpolator,color='#f98f60', alpha=0.3, label='Confidence Interval (95%)')
    # Add labels and title
    ax.set(xlabel=f'{p_name} ({p_units})', ylabel='NPV15 (MM USD)', xlim=(lower,upper))
    ax.grid(visible=True)
    ax.ticklabel_format(axis='x', style='sci')
    # Show the legend
    fig.tight_layout()
    ax.legend(loc='upper left')
    # Save the plot
    fig.savefig(os.path.join(figures_path, f'Target_{p_name}.tiff'), dpi=600)
