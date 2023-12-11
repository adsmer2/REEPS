import numpy as np 
import scipy as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import qsdsan as qs
import biosteam as bst
import os

from systems import create_system
from model import *

def analysis_uncertainty(sys, fununit, num_samples, figures_path):
    # Create model
    # -------------------------
    model_uncertainty = create_model(sys, fununit, 'all', 'no') # set system, functional unit, and type of parameter for analysis

    np.random.seed(3221) # setting the seed ensures you get the same sample

    samples = model_uncertainty.sample(N=num_samples, rule='L')
    model_uncertainty.load_samples(samples)
    model_uncertainty.evaluate()

    # 2-D kernel density plot with box/whisker in margins (for showing uncertainty for two indicators)
    # --------------------------
    for i in range(0, len(model_uncertainty.metrics)):
        if 'GWP100' in model_uncertainty.metrics[i].name:
            GWP_result = model_uncertainty.metrics[i]
        elif 'NPV15' in model_uncertainty.metrics[i].name:
            NPV_result = model_uncertainty.metrics[i]
        elif 'MSP' in model_uncertainty.metrics[i].name:
            MSP_result = model_uncertainty.metrics[i]

    xdata = model_uncertainty.table.loc[:,('TEA', 'NPV15 [MM USD/kg]')][lambda x:x>-4000] 
    ydata = model_uncertainty.table.loc[:,('LCA', 'Global warming [kg CO2-Eq/kg PG]')]
    plt.style.use('default')
    fig = sns.JointGrid()
    sns.kdeplot(x = xdata, y = ydata, fill=True, color= '#79bf82', ax = fig.ax_joint)
    sns.boxplot(x = xdata, ax = fig.ax_marg_x, color= '#79bf82')
    sns.boxplot(y = ydata, ax = fig.ax_marg_y, color= '#79bf82')
    fig.set_axis_labels(xlabel=f'NPV15 (MM USD)', ylabel=f'Global Warming (kg CO\u2082-eq/kg {fununit})')
    # fig, ax = qs.stats.plot_uncertainties(model_uncertainty, x_axis=NPV_result, y_axis=GWP_result, kind='kde-box', center_kws={'fill': True, 'color': '#79bf82'}, margin_kws={'color': '#79bf82'}) 
    # ax0, ax1, ax2 = fig.axes # KDE, top box, right box
    # ax0.set(xlabel=f'NPV (MM USD)', ylabel=f'Global Warming (kg CO2 Eq./kg {fununit})')
    fig.savefig(os.path.join(figures_path, f'kde_uncertainty_{fununit}.tiff'), dpi=600)
    
    return model_uncertainty


# # Data for five parameters
# # Data for five parameters
# param1_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,0]
# param2_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,1]
# param3_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,2]
# param4_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,3]
# param5_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,4]

# # Baseline values for each parameter
# baseline_values_all = model_uncertainty.metrics_at_baseline()['LCA',:]
# baseline_values = [baseline_values_all.iloc[0], 
# baseline_values_all.iloc[1], 
# baseline_values_all.iloc[2], 
# baseline_values_all.iloc[3], 
# baseline_values_all.iloc[4]]

# # Create subplots with shared y-axis
# # Create figures the correct size for publication
# aspect_ratio_LtoW = 6/4
# cm_to_in = 1/2.54  # centimeters in inches
# width_one_col = 8.3 # cm. Width for a one column figure
# width_two_col = 17.1 # cm. Width for a two column figure
# max_length = 23.3 # cm. The maximum lenght a figure can be
# plt.style.use('default')
# fig, axs = plt.subplots(5, sharey=False, figsize=(width_two_col*cm_to_in, width_two_col*aspect_ratio_LtoW*cm_to_in))

# # Extend ScalarFormatter
# class MyScalarFormatter(ScalarFormatter):
#     # Override '_set_format' with your own
#     def _set_format(self):
#         self.format = '%.1f'  # Show 1 decimal

# # red f1777f
# # blue 60c1cf
# # green 79bf82
# # orange f98f60
# # purple a280b9
# # gray 90918e
# # yellow f3c354
# # black 403a48
# colors = ['#f1777f', '#60c1cf', '#79bf82', '#f98f60', '#a280b9']

# # KDE plots for each parameter
# for i, (param_values, baseline, color) in enumerate(zip([param1_values, param2_values, param3_values, param4_values, param5_values], baseline_values, colors)):
#     sns.kdeplot(param_values, ax=axs[i], label=f'Parameter {i+1}', fill=True, color=color)

#     # Calculate the midpoint of the y-axis range
#     y_midpoint = 0.5 * (axs[i].get_ylim()[0] + axs[i].get_ylim()[1])
#     axs[i].scatter(baseline, y_midpoint, marker='o', label='Baseline', color=color)  # Adjusted y-coordinate
    
#     # Use Custom formatter for set number of decimals and scientific notation
#     custom_formatter = MyScalarFormatter(useMathText=True)
#     axs[i].xaxis.set_major_formatter(custom_formatter)
#     axs[i].xaxis.major.formatter.set_powerlimits((0,0))

#     # Remove left, right, and top spines for all subplots
#     axs[i].spines['left'].set_color('none')
#     axs[i].spines['right'].set_color('none')
#     axs[i].spines['top'].set_color('none')
#     axs[i].yaxis.tick_left()
#     axs[i].yaxis.set_ticks_position('none')  # Remove y-axis ticks
#     axs[i].set_yticklabels([])  # Remove y-axis labels
#     axs[i].set_ylabel('')  # Remove y-axis label

# # Apply tight layout for better spacing
# fig.tight_layout()

# # Show the plot
# plt.show()
