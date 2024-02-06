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

# def analysis_uncertainty(sys, fununit, num_samples, figures_path):
#     # Create model
#     # -------------------------
#     model_uncertainty = create_model(sys, fununit, 'all', 'no') # set system, functional unit, and type of parameter for analysis

#     np.random.seed(3221) # setting the seed ensures you get the same sample

#     samples = model_uncertainty.sample(N=num_samples, rule='L')
#     model_uncertainty.load_samples(samples)
#     model_uncertainty.evaluate()

#     # Creates figure of all 21 parameters
#     # Data for 21 parameters
#     paramLCA_values = model_uncertainty.table.loc[:, 'LCA']
#     paramTEA_values_all = model_uncertainty.table.loc[:, 'TEA']
#     paramTEA_values = paramTEA_values_all.loc[:, ('NPV15 [MM USD]', 'IRR [%]', 'MSP [USD/kg REO]')]
#     param_combined = pd.concat([paramLCA_values, paramTEA_values], axis=1)
#     series_list = [param_combined[:][col] for col in param_combined[:].columns]

#     # Baseline values for each parameter
#     baseline_values_LCA = model_uncertainty.metrics_at_baseline()['LCA', :]
#     baseline_values_TEA = model_uncertainty.metrics_at_baseline()['TEA', :]

#     # Sample parameter strings with names and units
#     param_strings_LCA = model_uncertainty.table.loc[:, 'LCA'].columns
#     param_strings_TEA_all = model_uncertainty.table.loc[:, 'TEA'].columns
#     param_strings_TEA = [i for i in param_strings_TEA_all if 'NPV15' in i or 'IRR' in i or 'MSP' in i]

#     # Extract parameter names and units
#     param_names_LCA = [param.split(' [')[0] for param in param_strings_LCA]
#     param_units_LCA = [param.split(' [')[1][:-1] for param in param_strings_LCA]  # Remove the trailing ']'
#     param_names_TEA = [param.split(' [')[0] for param in param_strings_TEA]
#     param_units_TEA = [param.split(' [')[1][:-1] for param in param_strings_TEA]  # Remove the trailing ']'

#     param_names = param_names_LCA + param_names_TEA
#     param_units = param_units_LCA + param_units_TEA

#     # Create figures the correct size for publication
#     aspect_ratio_LtoW = 2.8 # largest can be 2.807
#     cm_to_in = 1/2.54  # centimeters in inches
#     width_one_col = 8.3 # cm. Width for a one column figure
#     width_two_col = 17.1 # cm. Width for a two column figure
#     max_length = 23.3 # cm. The maximum lenght a figure can be
#     # Create subplots with shared y-axis
#     plt.style.use('default')
#     font = {'family': 'Calibri', 'size': 10}
#     plt.rc('font', **font)
#     fig, axs = plt.subplots(7, 3, figsize=(width_two_col*cm_to_in, 15))

#     # Extend ScalarFormatter
#     class MyScalarFormatter(ScalarFormatter):
#         # Override '_set_format' with your own
#         def _set_format(self):
#             self.format = '%.2f'  # Show 1 decimal

#     colors = sns.color_palette(palette='husl', n_colors=21).as_hex()[:]

#     # Flatten the axs array to iterate over subplots
#     axs_flat = axs.flatten()

#     # Loop over each parameter and create a KDE plot in the corresponding subplot
#     for i, (param_values, baseline, color) in enumerate(zip(series_list, pd.concat([baseline_values_LCA, baseline_values_TEA]), colors)):
#         row = i // 3
#         col = i % 3

#         sns.kdeplot(param_values, ax=axs[row, col], label=f'Parameter {i+1}', fill=True, color=color)

#         # Calculate the midpoint of the y-axis range
#         y_midpoint = 0.5 * (axs[row, col].get_ylim()[0] + axs[row, col].get_ylim()[1])

#         # Use the calculated y-coordinate for the vertical line with distinct color
#         # axs[row, col].axvline(baseline, ymin=0, ymax=0.2, linestyle='-', color=color)
#         axs[row, col].scatter(baseline, y_midpoint, marker='o', label='Baseline', color=color, s = 10)

#         # Set x-axis label to the right side with two lines
#         axs[row, col].set_xlabel(f'{param_names[i]}\n({param_units[i]})', ha='left', color=color)
#         axs[row, col].xaxis.set_label_coords(0.01, -0.35)

#         # Use Custom formatter for set number of decimals and scientific notation
#         custom_formatter = MyScalarFormatter(useMathText=True)
#         axs[row, col].xaxis.set_major_formatter(custom_formatter)
#         axs[row, col].xaxis.major.formatter.set_powerlimits((0, 0))

#         # Remove left, right, and top spines for all subplots
#         axs[row, col].spines['left'].set_color('none')
#         axs[row, col].spines['right'].set_color('none')
#         axs[row, col].spines['top'].set_color('none')
#         axs[row, col].yaxis.tick_left()
#         axs[row, col].yaxis.set_ticks_position('none')  # Remove y-axis ticks
#         axs[row, col].set_yticklabels([])  # Remove y-axis labels
#         axs[row, col].set_ylabel('')  # Remove y-axis label

#     # Apply tight layout for better spacing
#     fig.tight_layout()

#     # Show the plot
#     fig.savefig(os.path.join(figures_path, f'kde_uncertainty_{fununit}.tiff'), dpi=600)


def analysis_uncertainty(sys, fununit, num_samples, figures_path):
    # Create model
    # -------------------------
    model_uncertainty = create_model(sys, fununit, 'all', 'no') # set system, functional unit, and type of parameter for analysis

    np.random.seed(3221) # setting the seed ensures you get the same sample

    samples = model_uncertainty.sample(N=num_samples, rule='L')
    model_uncertainty.load_samples(samples)
    model_uncertainty.evaluate()
    # Creates figure of key 10 parameters for paper
    # ------------------------------------
    # Data for five parameters
    param1_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,0]
    param2_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,1]
    param3_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,2]
    param4_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,4]
    param5_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,6]
    param6_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,8]
    param7_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,9]
    param8_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,10]
    param9_values = model_uncertainty.table.loc[:,'LCA'].iloc[:,11]
    param10_values = model_uncertainty.table.loc[:,'TEA'].loc[:,'NPV15 [MM USD]']

    # Baseline values for each parameter
    baseline_values_all = model_uncertainty.metrics_at_baseline()['LCA',:]
    baseline_values = [baseline_values_all.iloc[0], 
    baseline_values_all.iloc[1], 
    baseline_values_all.iloc[2], 
    baseline_values_all.iloc[4], 
    baseline_values_all.iloc[6],
    baseline_values_all.iloc[8],
    baseline_values_all.iloc[9],
    baseline_values_all.iloc[10],
    baseline_values_all.iloc[11],
    model_uncertainty.metrics_at_baseline()['TEA',:].loc['NPV15 [MM USD]']
    ]

    # Sample parameter strings with names and units
    param_strings = [model_uncertainty.table.loc[:,'LCA'].iloc[:,0].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,1].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,2].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,4].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,6].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,8].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,9].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,10].name,
        model_uncertainty.table.loc[:,'LCA'].iloc[:,11].name,
        model_uncertainty.table.loc[:,'TEA'].loc[:,'NPV15 [MM USD]'].name
    ]

    # Extract parameter names and units
    param_names = [param.split(' [')[0] for param in param_strings]
    param_units = [param.split(' [')[1][:-1] for param in param_strings]  # Remove the trailing ']'

    # Create subplots with shared y-axis
    # Create figures the correct size for publication
    aspect_ratio_LtoW = 2.8 # largest can be 2.807
    cm_to_in = 1/2.54  # centimeters in inches
    width_one_col = 8.3 # cm. Width for a one column figure
    width_two_col = 17.1 # cm. Width for a two column figure
    max_length = 23.3 # cm. The maximum lenght a figure can be
    plt.style.use('default')
    # Manually change the text font and size
    font = {'family': 'Calibri', 'size': 8}
    plt.rc('font', **font)
    fig, axs = plt.subplots(5, 2, figsize=(width_one_col/2.54, 6.5))

    # Extend ScalarFormatter
    class MyScalarFormatter(ScalarFormatter):
        # Override '_set_format' with your own
        def _set_format(self):
            self.format = '%.1f'  # Show 1 decimal

    colors = sns.color_palette(palette='husl', n_colors=10).as_hex()[:]

    # Flatten the axs array to iterate over subplots
    axs_flat = axs.flatten()

    # KDE plots for each parameter
    for i, (param_values, baseline, color) in enumerate(zip([param1_values, param2_values, param3_values, param4_values, param5_values, param6_values, param7_values, param8_values, param9_values, param10_values], baseline_values, colors)):
        row = i // 2
        col = i % 2

        sns.kdeplot(param_values, ax=axs[row, col], label=f'Parameter {i+1}', fill=True, color=color)

        # Calculate the midpoint of the y-axis range
        y_midpoint = 0.5 * (axs[row, col].get_ylim()[0] + axs[row, col].get_ylim()[1])
        axs[row, col].scatter(baseline, y_midpoint, marker='o', label='Baseline', color=color, s = 10)  # Adjusted y-coordinate

        # Use the calculated y-coordinate for the vertical line with distinct color
        # axs[row, col].axvline(baseline, ymin=0, ymax=0.2, linestyle='-', color=color)

        # Set x-axis label to the right side with two lines
        axs[row, col].set_xlabel(f'{param_names[i]}\n({param_units[i]})', ha='left', color=color)
        axs[row, col].xaxis.set_label_coords(-0.01, -0.42)

        # Use Custom formatter for set number of decimals and scientific notation
        custom_formatter = MyScalarFormatter(useMathText=True)
        axs[row, col].xaxis.set_major_formatter(custom_formatter)
        axs[row, col].xaxis.major.formatter.set_powerlimits((0, 0))

        # Remove left, right, and top spines for all subplots
        axs[row, col].spines['left'].set_color('none')
        axs[row, col].spines['right'].set_color('none')
        axs[row, col].spines['top'].set_color('none')
        axs[row, col].yaxis.tick_left()
        axs[row, col].yaxis.set_ticks_position('none')  # Remove y-axis ticks
        axs[row, col].set_yticklabels([])  # Remove y-axis labels
        axs[row, col].set_ylabel('')  # Remove y-axis label

    # Apply tight layout for better spacing
    fig.tight_layout()

    # Show the plot
    fig.savefig(os.path.join(figures_path, f'kde_uncertainty_paper_{fununit}.tiff'), dpi=600)
    return model_uncertainty
    



# # 2-D kernel density plot with box/whisker in margins (for showing uncertainty for two indicators)
# # --------------------------
# for i in range(0, len(model_uncertainty.metrics)):
#     if 'GWP100' in model_uncertainty.metrics[i].name:
#         GWP_result = model_uncertainty.metrics[i]
#     elif 'NPV15' in model_uncertainty.metrics[i].name:
#         NPV_result = model_uncertainty.metrics[i]
#     elif 'MSP' in model_uncertainty.metrics[i].name:
#         MSP_result = model_uncertainty.metrics[i]

# xdata = model_uncertainty.table.loc[:,('TEA', 'NPV15 [MM USD/kg]')][lambda x:x>-4000] 
# ydata = model_uncertainty.table.loc[:,('LCA', 'Global warming [kg CO2-Eq/kg PG]')]
# plt.style.use('default')
# fig = sns.JointGrid()
# sns.kdeplot(x = xdata, y = ydata, fill=True, color= '#79bf82', ax = fig.ax_joint)
# sns.boxplot(x = xdata, ax = fig.ax_marg_x, color= '#79bf82')
# sns.boxplot(y = ydata, ax = fig.ax_marg_y, color= '#79bf82')
# fig.set_axis_labels(xlabel=f'NPV15 (MM USD)', ylabel=f'Global Warming (kg CO\u2082-eq/kg {fununit})')
# # fig, ax = qs.stats.plot_uncertainties(model_uncertainty, x_axis=NPV_result, y_axis=GWP_result, kind='kde-box', center_kws={'fill': True, 'color': '#79bf82'}, margin_kws={'color': '#79bf82'}) 
# # ax0, ax1, ax2 = fig.axes # KDE, top box, right box
# # ax0.set(xlabel=f'NPV (MM USD)', ylabel=f'Global Warming (kg CO2 Eq./kg {fununit})')
# fig.savefig(os.path.join(figures_path, f'2param_uncertainty_{fununit}.tiff'), dpi=600)

# return model_uncertainty