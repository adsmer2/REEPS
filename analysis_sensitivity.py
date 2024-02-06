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

def analysis_sensitivity(sys, fununit, parameter, num_samples, figures_path):
    # Sensitivity Analysis
    # -------------------------
    # use above model to get correlations 
    model_sensitivity = create_model(sys, fununit, parameter, 'no') # set system, functional unit, and type of parameter for analysis

    np.random.seed(3221) # setting the seed ensures you get the same sample

    samples = model_sensitivity.sample(N=num_samples, rule='L')
    model_sensitivity.load_samples(samples)
    model_sensitivity.evaluate()
    r_df, p_df = qs.stats.get_correlations(model_sensitivity, kind='Spearman')
    
    # sort the parameters by alphabetically by unit ID then parameter name
    lst1 = [i.element for i in model_sensitivity.get_parameters()]
    lst2 = [i.name for i in model_sensitivity.get_parameters()]
    key_parameters = sorted(list(zip(lst1, lst2)), key=lambda x: x[0])

    # Manually make changes to bubble plot built off functions from qsdsan.stats
    def _update_df_names(df, columns=True, index=True):
        new_df = df.copy()

        if columns and not new_df.columns.empty:
            try:
                iter(new_df.columns)
                new_df.columns = [i[-1].split(' [')[0] for i in new_df.columns]
            except: pass

        if index and not new_df.index.empty:
            try:
                iter(new_df.index)
                new_df.index = [i[-1].split(' [')[0] for i in new_df.index]
            except: pass

        return new_df

    def _update_input(input_val, default_val):
        if input_val is None:
            return default_val
        else:
            try:
                iter(input_val)
                if len(input_val)==0: # empty iterable
                    return default_val
                return input_val if not isinstance(input_val, str) else (input_val,)
            except:
                return (input_val,)

    df = _update_df_names(r_df)

    param_names = _update_input(np.array(key_parameters)[:,1], df.index)
    param_names = param_names if isinstance(param_names[0], str) \
                                else [p.name for p in param_names]
    metric_names = _update_input(model_sensitivity.metrics, df.columns)
    metric_names = metric_names if isinstance(metric_names[0], str) \
                                else [m.name for m in metric_names]

    df = df[metric_names].loc[param_names]

    corr_df = df.stack(dropna=False).reset_index()
    corr_df.rename(columns={'level_0': 'parameter', 'level_1': 'metric',
                            0: 'Sign'}, inplace=True)
    corr_df['Correlation'] = corr_df['Sign'].abs()

    # correlation dataframe for heatmap plot
    corr_df2 = df.stack(dropna=False).reset_index()
    corr_df2.rename(columns={'level_0': 'parameter', 'level_1': 'metric',
                            0: 'Sign'}, inplace=True)
    corr_df2['Correlation'] = corr_df['Sign']

    # # Fix the sign of the NLTP CFs
    # condition = (corr_df2['metric'] == 'Natural land transformation')
    # corr_df2.loc[condition, 'Correlation'] = -corr_df2.loc[condition, 'Correlation']

    # Remove all metrics that aren't economic because enviornmental impacts are unchanged by price changes
    if parameter == 'contextual':
        corr_df2 = corr_df2[~corr_df2['metric'].isin(metric_names[3:])] 
        pivot_corr_df = corr_df2.pivot(index="parameter", columns="metric", values="Correlation")
        pivot_corr_df = pivot_corr_df.reindex(columns=['MSP',
        'NPV15',
        'IRR',
        ])
        columns = ['Minimum Selling Price',
        'Net Present Value',
        'Internal Rate of Return'
        ]
        pivot_corr_df.columns = columns

    if parameter == 'technological':
        # List of exceptions
        exceptions = ['NPV15', 'IRR', 'MSP']
        # Custom function to apply title() only if not in exceptions
        def custom_title(string):
            return string.title() if string not in exceptions else string
        # Apply the custom_title function to the desired column
        corr_df2['metric'] = corr_df2['metric'].apply(custom_title)

        pivot_corr_df = corr_df2.pivot(index="parameter", columns="metric", values="Correlation")
        pivot_corr_df = pivot_corr_df.sort_index(key=lambda x: [i.split('(')[1] for i in x])
        # Need to put the categories in the correct order for plotting
        pivot_corr_df = pivot_corr_df.reindex(columns=['MSP',
        'Acidification Terrestrial',
        'Climate Change',
        'Ecotoxicity Freshwater',
        'Ecotoxicity Marine',
        'Ecotoxicity Terrestrial',
        'Energy Resources',
        'Eutroph. Freshwater',
        'Eutroph. Marine',
        'Human Toxicity Carc.',
        'Human Toxicity N-Carc.', # watch out for capitalization after punctuation
        'Ionising Radiation',
        'Land Use',
        'Meterial Resources',
        'Ozone Depletion',
        'Particulate Matter',
        'Photochemical Ox. Human Health',
        'Photochemical Ox. Ecosystems',
        'Water Use',
        'NPV15',
        'IRR'
        ])
        pivot_corr_df.drop(['NPV15','IRR'], inplace=True, axis=1)
        columns = ['Minimum Selling Price',
        'Acidification Terrestrial',
        'Climate Change',
        'Ecotoxicity Freshwater',
        'Ecotoxicity Marine',
        'Ecotoxicity Terrestrial',
        'Energy Resources',
        'Eutroph. Freshwater',
        'Eutroph. Marine',
        'Human Toxicity Carc.',
        'Human Toxicity N-Carc.', # watch out for capitalization after punctuation
        'Ionising Radiation',
        'Land Use',
        'Meterial Resources',
        'Ozone Depletion',
        'Particulate Matter',
        'Photochemical Ox. Human Health',
        'Photochemical Ox. Ecosystems',
        'Water Use'
        ]
        pivot_corr_df.columns = columns
        

    # Manually change the text font and size
    plt.style.use('default')
    font = {'family': 'Calibri', 'size': 8}
    plt.rc('font', **font)

    # if parameter == 'technological':
    #     fig, ax = plt.subplots(figsize=(6.73, 6.73*0.8))
    #     sns.heatmap(data=pivot_corr_df.T, vmin = -1, vmax=1, cmap=sns.diverging_palette(356, 200, s=80, l=50, as_cmap=True))
    # if parameter == 'contextual':
    #     fig, ax = plt.subplots(figsize=(6.73*0.82, 6.73*0.4))
    #     sns.heatmap(data=pivot_corr_df.T, vmin = -1, vmax=1, cmap=sns.diverging_palette(356, 200, s=80, l=50, as_cmap=True))

    if parameter == 'technological':
        fig, ax = plt.subplots(figsize=(6.73, 6.73*0.8))
        sns.heatmap(data=pivot_corr_df.T, vmin = -1, vmax=1, cmap=sns.diverging_palette(220, 356, s=90, l=60, as_cmap=True))
    if parameter == 'contextual':
        fig, ax = plt.subplots(figsize=(6.73*0.95, 6.73*0.4))
        sns.heatmap(data=pivot_corr_df.T, vmin = -1, vmax=1, cmap=sns.diverging_palette(220, 356, s=90, l=60, as_cmap=True))

    ax.set(xlabel='', ylabel='')
    ax.tick_params(axis='y', labelrotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(figures_path, f'Sensitivity_heatmap_{parameter}_{fununit}_FULL.tiff'), dpi=600)

    return model_sensitivity


# Morris OAT Analysis
# -------------------------
# model_morris = create_model(sys, fununit, parameter, target) # set system, functional unit, and type of parameter for analysis

# inputs = qs.stats.define_inputs(model_morris)
# samples_morris = qs.stats.generate_samples(inputs, kind='Morris', N=10, seed=554) # num_levels=num_levels. Default is 4

# model_morris.load_samples(samples_morris)
# model_morris.evaluate()

# dct = qs.stats.morris_analysis(model_morris, inputs, seed=554, nan_policy='fill_mean', file=os.path.join(results_path, f'Morris_Sensitivity_{fununit}_{parameter}.xlsx'))
# fig, ax = qs.stats.plot_morris_results(dct, metric=model.metrics[0]) 
# fig.savefig(os.path.join(figures_path, f'Morris_Sensitivity_{fununit}_{parameter}.png'), dpi=300)


# # Bubble sensitivity plot
# # --------------------------
# # use above model to get correlations 
# model_sensitivity = create_model(sys, fununit, parameter, 'no') # set system, functional unit, and type of parameter for analysis

# np.random.seed(3221) # setting the seed ensures you get the same sample

# samples = model_sensitivity.sample(N=num_samples, rule='L')
# model_sensitivity.load_samples(samples)
# model_sensitivity.evaluate()
# r_df, p_df = qs.stats.get_correlations(model_sensitivity, kind='Spearman')

# # sort the parameters by alphabetically by unit ID then parameter name
# lst1 = [i.element for i in model_sensitivity.get_parameters()]
# lst2 = [i.name for i in model_sensitivity.get_parameters()]
# key_parameters = sorted(list(zip(lst1, lst2)), key=lambda x: x[0])

# # # Filter out parameters that only meet a certain threshold
# # def filter_parameters(model, df, threshold):
# #     new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
# #     filtered = new_df.dropna(how='all')
# #     param_dct = {p.name_with_units:p for p in model.get_parameters()}
# #     parameters = set(param_dct[i[1]] for i in filtered.index)
# #     return list(parameters)
# # key_parameters2 = filter_parameters(model_sensitivity, r_df, threshold=0.15) # Only want parameters with Spearman's rho >= 0.4 or <= -0.4

# # Manually make changes to bubble plot built off functions from qsdsan.stats
# def _update_df_names(df, columns=True, index=True):
#     new_df = df.copy()

#     if columns and not new_df.columns.empty:
#         try:
#             iter(new_df.columns)
#             new_df.columns = [i[-1].split(' [')[0] for i in new_df.columns]
#         except: pass

#     if index and not new_df.index.empty:
#         try:
#             iter(new_df.index)
#             new_df.index = [i[-1].split(' [')[0] for i in new_df.index]
#         except: pass

#     return new_df

# def _update_input(input_val, default_val):
#     if input_val is None:
#         return default_val
#     else:
#         try:
#             iter(input_val)
#             if len(input_val)==0: # empty iterable
#                 return default_val
#             return input_val if not isinstance(input_val, str) else (input_val,)
#         except:
#             return (input_val,)

# df = _update_df_names(r_df)

# filtered_unit_name = [i.element for i in key_parameters2]
# filtered_param_name = [i.name for i in key_parameters2]
# key_parameters2 = sorted(list(zip(filtered_unit_name, filtered_param_name)), key=lambda x: x[0])
# param_names = _update_input(np.array(key_parameters2)[:,1], df.index)
# # param_names = [i.name for i in key_parameters2]
# param_names = param_names if isinstance(param_names[0], str) \
#                             else [p.name for p in param_names]
# metric_names = _update_input(model_sensitivity.metrics, df.columns)
# metric_names = metric_names if isinstance(metric_names[0], str) \
#                             else [m.name for m in metric_names]

# df = df[metric_names].loc[param_names]


# corr_df = df.stack(dropna=False).reset_index()
# corr_df.rename(columns={'level_0': 'parameter', 'level_1': 'metric',
#                         0: 'Sign'}, inplace=True)
# corr_df['Correlation'] = corr_df['Sign'].abs()

# # correlation dataframe for heatmap plot
# corr_df2 = df.stack(dropna=False).reset_index()
# corr_df2.rename(columns={'level_0': 'parameter', 'level_1': 'metric',
#                         0: 'Sign'}, inplace=True)
# corr_df2['Correlation'] = corr_df['Sign']

# # make DataFrame
# data = {'Sign': corr_df['Sign']}
# df_gpt = pd.DataFrame(data)

# # Function to categorize values
# def categorize_sign(value):
#     if value < 0:
#         return '$-$'
#     elif value > 0:
#         return '$+$'
#     else:
#         return 'Zero'

# # Apply the function to create a new column
# corr_df['Sign'] = df_gpt['Sign'].apply(categorize_sign)

# # Remove all metrics that aren't economic because enviornmental impacts are unchanged by price changes
# if parameter == 'contextual':
#     corr_df = corr_df[~corr_df['metric'].isin(metric_names[3:])] 

# # Begin plotting
# def _plot_corr_bubble(corr_df, ratio, **kwargs):
#     plt.style.use('default')

#     if parameter == 'technological':
#         margin_x = kwargs['margin_x'] if 'margin_x' in kwargs.keys() else 0.05
#         margin_y = kwargs['margin_y'] if 'margin_y' in kwargs.keys() else 0.05
#         kwargs = {i: kwargs[i] for i in kwargs.keys() if 'margin' not in i}

#         keys = ('height', 'aspect', 'sizes', 'size_norm', 'edgecolor') # , 'palette'
#         values = (9+ratio, 1, (0, 1000), (0, 2.5), '0.5') # , ['#60c1cf', '#f1777f']
#     elif parameter == 'contextual':
#         margin_x = kwargs['margin_x'] if 'margin_x' in kwargs.keys() else 0.2/ratio
#         margin_y = kwargs['margin_y'] if 'margin_y' in kwargs.keys() else 0.05
#         kwargs = {i: kwargs[i] for i in kwargs.keys() if 'margin' not in i}

#         keys = ('height', 'aspect', 'sizes', 'size_norm', 'edgecolor') # , 'palette'
#         values = (7+ratio, 0.7, (0, 1000), (0, 2.5), '0.5') # , ['#60c1cf', '#f1777f']
#     else: 
#         RuntimeError(f'parameter={parameter} is not "technological" or "contextual". Please define as one of these two.')

#     for num, k in enumerate(keys):
#         kwargs.setdefault(keys[num], values[num])

#     g = sns.relplot(data=corr_df, x='parameter', y='metric',
#                     hue='Sign', size='Correlation', palette= {'$+$':'#60c1cf', '$-$':'#f1777f'}, **kwargs)

#     g.set(xlabel='', ylabel='', aspect=1)
#     g.ax.margins(x=margin_x, y=margin_y)

#     for label in g.ax.get_xticklabels():
#         label.set_rotation(90)

#     # for artist in g.legend.legendHandles:
#     #     artist.set_edgecolor('1')

#     for key in g.ax.spines.keys():
#         g.ax.spines[key].set(color='k', linewidth=0.5, visible=True)

#     g.ax.grid(True, which='major', color='k',linestyle='--', linewidth=0.3)
#     g.tight_layout()
    
#     if parameter == 'technological':
#         sns.move_legend(g, 'center right', bbox_to_anchor=(0.65, 0.55))
#     elif parameter == 'contextual':
#         sns.move_legend(g, 'center right', bbox_to_anchor=(0.85, 0.55))
#     else: 
#         RuntimeError(f'parameter={parameter} is not "technological" or "contextual". Please define as one of these two.')

#     return g

# g = _plot_corr_bubble(corr_df, len(metric_names)/len(param_names))
# g.savefig(os.path.join(figures_path, f'bubble_sensitivity_{fununit}_{parameter}.tiff'), dpi=600)

# g = _plot_corr_bubble(corr_df, len(metric_names)/len(param_names))
# g.savefig(os.path.join(figures_path, f'bubble_sensitivity_{fununit}_{parameter}.tiff'), dpi=600)



# # the most recent version with figures that fit the paper and fixes
# # Sensitivity Analysis
# # -------------------------
# # use above model to get correlations 
# model_sensitivity = create_model(sys, fununit, parameter, 'no') # set system, functional unit, and type of parameter for analysis

# np.random.seed(3221) # setting the seed ensures you get the same sample

# samples = model_sensitivity.sample(N=num_samples, rule='L')
# model_sensitivity.load_samples(samples)
# model_sensitivity.evaluate()
# r_df, p_df = qs.stats.get_correlations(model_sensitivity, kind='Spearman')


# # # sort the parameters by alphabetically by unit ID then parameter name
# # lst1 = [i.element for i in model_sensitivity.get_parameters()]
# # lst2 = [i.name for i in model_sensitivity.get_parameters()]
# # key_parameters = sorted(list(zip(lst1, lst2)), key=lambda x: x[0])

# # Filter out parameters that only meet a certain threshold
# def filter_parameters(model, df, threshold):
#     new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
#     filtered = new_df.dropna(how='all')
#     param_dct = {p.name_with_units:p for p in model.get_parameters()}
#     parameters = set(param_dct[i[1]] for i in filtered.index)
#     return list(parameters)
# key_parameters2 = filter_parameters(model_sensitivity, r_df, threshold=0.15) # Only want parameters with Spearman's rho >= 0.4 or <= -0.4

# # Manually make changes to bubble plot built off functions from qsdsan.stats
# def _update_df_names(df, columns=True, index=True):
#     new_df = df.copy()

#     if columns and not new_df.columns.empty:
#         try:
#             iter(new_df.columns)
#             new_df.columns = [i[-1].split(' [')[0] for i in new_df.columns]
#         except: pass

#     if index and not new_df.index.empty:
#         try:
#             iter(new_df.index)
#             new_df.index = [i[-1].split(' [')[0] for i in new_df.index]
#         except: pass

#     return new_df

# def _update_input(input_val, default_val):
#     if input_val is None:
#         return default_val
#     else:
#         try:
#             iter(input_val)
#             if len(input_val)==0: # empty iterable
#                 return default_val
#             return input_val if not isinstance(input_val, str) else (input_val,)
#         except:
#             return (input_val,)

# df = _update_df_names(r_df)

# filtered_unit_name = [i.element for i in key_parameters2]
# filtered_param_name = [i.name for i in key_parameters2]
# key_parameters2 = sorted(list(zip(filtered_unit_name, filtered_param_name)), key=lambda x: x[0])
# param_names = _update_input(np.array(key_parameters2)[:,1], df.index)
# # param_names = [i.name for i in key_parameters2]
# param_names = param_names if isinstance(param_names[0], str) \
#                             else [p.name for p in param_names]
# metric_names = _update_input(model_sensitivity.metrics, df.columns)
# metric_names = metric_names if isinstance(metric_names[0], str) \
#                             else [m.name for m in metric_names]

# df = df[metric_names].loc[param_names]


# corr_df = df.stack(dropna=False).reset_index()
# corr_df.rename(columns={'level_0': 'parameter', 'level_1': 'metric',
#                         0: 'Sign'}, inplace=True)
# corr_df['Correlation'] = corr_df['Sign'].abs()

# # correlation dataframe for heatmap plot
# corr_df2 = df.stack(dropna=False).reset_index()
# corr_df2.rename(columns={'level_0': 'parameter', 'level_1': 'metric',
#                         0: 'Sign'}, inplace=True)
# corr_df2['Correlation'] = corr_df['Sign']

# # Fix the sign of the NLTP CFs
# condition = (corr_df2['metric'] == 'Natural land transformation')
# corr_df2.loc[condition, 'Correlation'] = -corr_df2.loc[condition, 'Correlation']

# # Remove all metrics that aren't economic because enviornmental impacts are unchanged by price changes
# if parameter == 'contextual':
#     corr_df2 = corr_df2[~corr_df2['metric'].isin(metric_names[3:])] 
#     pivot_corr_df = corr_df2.pivot(index="parameter", columns="metric", values="Correlation")

# if parameter == 'technological':
#     # List of exceptions
#     exceptions = ['NPV15', 'IRR', 'MSP']
#     # Custom function to apply title() only if not in exceptions
#     def custom_title(string):
#         return string.title() if string not in exceptions else string
#     # Apply the custom_title function to the desired column
#     corr_df2['metric'] = corr_df2['metric'].apply(custom_title)

#     pivot_corr_df = corr_df2.pivot(index="parameter", columns="metric", values="Correlation")
#     pivot_corr_df = pivot_corr_df.sort_index(key=lambda x: [i.split('(')[1] for i in x])
#     pivot_corr_df = pivot_corr_df.reindex(columns=['Ag. Land Occupation',
#     'Fossil Depletion',
#     'Freshwater Ecotoxicity',
#     'Freshwater Eut.',
#     'Global Warming',
#     'Human Toxicity',
#     'Ionising Radiation',
#     'Marine Ecotoxicity',
#     'Marine Eutrophication',
#     'Metal Depletion',
#     'Natural Land Transformation',
#     'Ozone Depletion',
#     'Particulate Matter Formation',
#     'Photochemical Oxidant Formation',
#     'Terrestrial Acidification',
#     'Terrestrial Ecotoxicity',
#     'Urban Land Occupation',
#     'Water Depletion',
#     'NPV15',
#     'IRR',
#     'MSP'])

# # Create figures the correct size for publication
# aspect_ratio_LtoW = 1.75
# cm_to_in = 1/2.54  # centimeters in inches
# width_one_col = 8.3 # cm. Width for a one column figure
# width_two_col = 17.1 # cm. Width for a two column figure
# max_length = 23.3 # cm. The maximum lenght a figure can be

# # Manually change the text font and size
# plt.style.use('default')
# font = {'family': 'Calibri', 'size': 8}
# plt.rc('font', **font)

# if parameter == 'technological':
#     fig, ax = plt.subplots(figsize=(11*cm_to_in, 6))
#     sns.heatmap(data=pivot_corr_df.T, vmin = -1, vmax=1, cmap=sns.diverging_palette(360, 220, s=100, l=50, as_cmap=True))
#     ax.set(xlabel='', ylabel='')
#     ax.tick_params(axis='y', labelrotation=0)
#     fig.tight_layout()
#     fig.savefig(os.path.join(figures_path, f'Sensitivity_heatmap_{parameter}_{fununit}.tiff'), dpi=600)
# if parameter == 'contextual':
#     fig, ax = plt.subplots(figsize=(6*cm_to_in, 2))
#     sns.heatmap(data=pivot_corr_df.T, vmin = -1, vmax=1, cmap=sns.diverging_palette(360, 220, s=100, l=50, as_cmap=True))
#     ax.set(xlabel='', ylabel='')
#     ax.tick_params(axis='y', labelrotation=0)
#     fig.tight_layout()
#     fig.savefig(os.path.join(figures_path, f'Sensitivity_heatmap_{parameter}_{fununit}.tiff'), dpi=600)