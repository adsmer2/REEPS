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



def analysis_MSP_contributions(sys, tea, fs_stream, fs_unit, fununit, lca_results, figures_path):
    """
    Calculations of contributors to MSP performed according to NREL methods


    important results tables:
    ----------------
    MSP_table                 -     high level MSP contributions
    
    MSP_unit_table            -     process section level level MSP contributions
    
    indicator_contributions   -     process section level indicator contributions
    """
    # Calculate the MSP for the system
    MSP = tea.solve_price(fs_stream.Ln2O3) # $/kg REO
    # Calculate Average Income Tax (NREL method)
    CF_table = tea.get_cashflow_table() # cashflow table (dataframe)
    discount_factor_CF = CF_table.iloc[:,14] # discount factor for each year 
    tax_CF = CF_table.iloc[:,10] # income tax ($)
    tax_CF_sum = np.sum(tax_CF*discount_factor_CF)*1000000 # discounted income tax over project lifespan ($)
    sales_CF = [sys.get_market_value(fs_stream.Ln2O3)]*len(tax_CF) # sales ($)
    sales_CF_sum = np.sum(sales_CF*discount_factor_CF) # total discounted sales over project lifespan ($)
    average_income_tax = tax_CF_sum/sales_CF_sum*MSP # NREL defined average income tax ($/kg REO)
    # Calculate other Contributors to MSP
    product_mflow = fs_stream.Ln2O3.F_mass*24*tea.operating_days # mass flow rate of REO product in kg/year
    gypsum_sales_MSP = -sys.get_market_value(fs_stream.gypsum)/product_mflow # $/kg REO. negative sign bc this is a byproduct not a cost.
    VOC_MSP = tea.VOC/product_mflow # variable operating cost ($/kg REO)
    FOC_MSP = tea.FOC/product_mflow # fixed operating cost ($/kg REO)
    capital_dep_MSP = tea.FCI/(tea.duration[1]-tea.duration[0])/product_mflow # capital depreciation ($/kg REO)
    # Calculate Average Return on Investment as defined by NREL
    average_ROI_MSP = MSP-np.sum([VOC_MSP, gypsum_sales_MSP, FOC_MSP, capital_dep_MSP, average_income_tax]) # $/kg REO
    # Get the high level process economics data to calculate contributions to MSP
    MSP_table_values = [MSP, VOC_MSP, gypsum_sales_MSP, FOC_MSP, capital_dep_MSP, average_income_tax, average_ROI_MSP]
    MSP_table_index = ['MSP', 'VOC', 'Gypsum Credit', 'FOC', 'Capital Depreciation', 'Average Income Tax', 'Average Return on Investment']
    MSP_table = pd.DataFrame(MSP_table_values)
    MSP_table.index = MSP_table_index
    MSP_table.columns = [''] # Manufacturing Costs (USD/kg REO)

    # Create a stacked bar chart
    # blue 82cfd0
    # blue/green 00a996
    # green 3ba459
    # yellow/green 8ead3e
    # dark green 007f3d
    # orange fcb813
    # brown 98876e
    # gray 403a48
    custom_colors = ['#82cfd0', '#403a48', '#fcb813', '#007f3d', '#8ead3e', '#3ba459'] # order of legend VOC, credit, FOC, cap deprec., tax, ROI

    # Create figures the correct size for publication
    aspect_ratio_LtoW = 1
    cm_to_in = 1/2.54  # centimeters in inches
    width_one_col = 8.3 # cm. Width for a one column figure
    width_two_col = 17.1 # cm. Width for a two column figure
    max_length = 23.3 # cm. The maximum lenght a figure can be

    plt.style.use('default')
    # Manually change the text font and size
    plt.style.use('default')
    font = {'family': 'Calibri', 'size': 8}
    plt.rc('font', **font)

    fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in/2, width_one_col*aspect_ratio_LtoW*cm_to_in))
    ax = MSP_table.drop('MSP', axis=0).T.plot(kind='bar', stacked=True, color=custom_colors, ax=ax)
    # Set labels and title
    ax.set_ylabel('Minimum Selling Price (USD/kg REO)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.6) ) # , ncol=len(combined_df.index)
    ax.axhline(y=0, color='black', linewidth=0.8)
    fig.tight_layout()
    # Show the plot
    fig.savefig(os.path.join(figures_path, f'MSP_Contributions.tiff'), dpi=600)

    # Gather process economics results for each unit and each inlet flow
    result = [] 
    for i in fs_unit:
        for j in i.ins:
            feed_cost = j.price*j.F_mass*24*tea.operating_days
            PandH = i.utility_cost*tea.operating_days*24 + sum([k for k in i.add_OPEX.values()])*tea.operating_days*24
            result.append([i.ID, j.price, j.F_mass, i.purchase_cost, feed_cost, PandH])
    econ_results_unit = pd.DataFrame(result) # process economics results by unit and by inlet flow
    econ_results_unit.index = econ_results_unit.iloc[:,0]
    econ_results_unit.columns = ['ID', 'Price', 'Mass Flow', 'Purchase Cost', 'Chemicals and Materials', 'Utilities']

    # Group by index and aggregate rows performing different math operations for each column
    econ_results_combined = econ_results_unit.groupby(econ_results_unit.index).agg({'Purchase Cost':'mean', 'Chemicals and Materials':'sum', 'Utilities':'mean'})
    econ_results_combined.loc['S1','Utilities'] = fs_unit.S1.utility_cost*tea.operating_days*24 # Do not want membrane cost here
    econ_results_combined.loc['S1','Chemicals and Materials'] = [i for i in fs_unit.S1.add_OPEX.values()][0]*24*tea.operating_days + econ_results_combined.loc['S1','Chemicals and Materials'] # Put membrane cost here
    econ_results_combined_2 = econ_results_combined.copy()

    # Create an index map to aggregate unit ID to process section
    index_mapping = {'U1': 'Leaching', 'F1': 'Leaching', 'F2': 'Leaching', 'M1': 'Leaching',
        'P1': 'Concentration', 'F3': 'Concentration', 'H1': 'Concentration',
        'P2': 'Refining', 'F4': 'Refining', 'H2': 'Refining', 'M2': 'Refining', 'M3': 'Refining',
        'S1': 'Separation', 'RS': 'Separation',
        'WT':'Wastewater Treatment', 'P3': 'Wastewater Treatment',
        'M4': 'REO Credit',
        'M0': 'PG Remediation Credit'}
    # Replace the index values based on the mapping
    econ_results_combined_2.index = econ_results_combined.index.to_series().replace(index_mapping)
    # Group by index and aggregate 'Value' using sum
    econ_results_grouped = econ_results_combined_2.groupby(econ_results_combined_2.index).agg({'Purchase Cost':'sum', 'Chemicals and Materials':'sum', 'Utilities':'sum'})
    econ_results_grouped.loc['Separation', 'Purchase Cost'] = fs_unit.S1.membrane_cost/4.28

    # Create the MSP_unit table
    MSP_unit_table = econ_results_grouped.copy()
    MSP_unit_table = MSP_unit_table.reindex(['Leaching', 'Concentration', 'Separation', 'Refining', 'Wastewater Treatment', 'REO Credit', 'PG Remediation Credit'])
    # insert a column for the gypsum credit
    gypsum_credit = [sys.get_market_value(fs_stream.gypsum), 0, 0, 0, 0, 0, 0]
    MSP_unit_table.insert(3, 'Gypsum Credit', gypsum_credit)
    # Convert econ_results_grouped into MSP contributions
    MSP_unit_table['Chemicals and Materials'] = MSP_unit_table['Chemicals and Materials'].multiply(1/product_mflow)
    MSP_unit_table['Utilities'] = MSP_unit_table['Utilities'].multiply(1/product_mflow)
    MSP_unit_table['Gypsum Credit'] = MSP_unit_table['Gypsum Credit'].multiply(-1/product_mflow)
    cap_recovery_charge = MSP_table.loc['Average Return on Investment'].values + MSP_table.loc['Average Income Tax'].values + MSP_table.loc['Capital Depreciation'].values
    MSP_unit_table['Fixed Operating Cost'] = MSP_unit_table.loc[:,'Purchase Cost'].multiply(4.28/tea.FCI*MSP_table.loc['FOC'].values[0])
    MSP_unit_table['Capital Recovery Charge'] = MSP_unit_table.loc[:,'Purchase Cost'].multiply(4.28/tea.FCI*cap_recovery_charge[0])
    MSP_unit_table.drop('Purchase Cost', inplace=True, axis=1)
    MSP_unit_table.drop('PG Remediation Credit', inplace=True, axis=0)

    # Create a stacked bar chart of unit level MSP contributors
    # blue 82cfd0
    # green 3ba459
    # orange fcb813
    # brown 98876e
    # gray 90918e
    # dark gray 403a48
    # dark blue 5973a6
    custom_colors = ['#82cfd0', '#5973a6', '#90918e', '#fcb813', '#3ba459']

    # Manually change the text font and size
    plt.style.use('default')
    font = {'family': 'Calibri', 'size': 11}
    plt.rc('font', **font)

    aspect_ratio_LtoW = 1.5 # Length/Width
    fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))
    ax = MSP_unit_table.drop('REO Credit', axis=0).plot(kind='bar', stacked=True, color=custom_colors, ax=ax)
    # Set labels and title
    ax.set_ylabel('Minimum Selling Price (USD/kg REO)')
    ax.set_xlabel('')
    ax.legend().remove()
    fig.legend(loc='upper left', bbox_to_anchor=(0.2, 1.01)) # , ncol=len(combined_df.index) ,bbox_to_anchor=(1.25, 0.65)
    # ax.legend(loc='lower right')
    ax.axhline(y=0, color='black', linewidth=0.8)
    fig.tight_layout()
    fig.subplots_adjust(top=0.75)
    # Show the plot
    fig.savefig(os.path.join(figures_path, f'MSP_Unit_Contributions.tiff'), dpi=600)


    # # Confirm TEA class methods match system class methods
    # tea.system.material_cost == feed_utility_table.loc[:,'Feed Cost'].sum() # Confirm Feed Cost
    # round(tea.utility_cost/1000000,6) == round(feed_utility_table.loc[:,'Utility Cost'].sum()/1000000,6) # Confirm Utility Cost
    # MSP_full_table.sum().sum() == MSP # Confirm overall MSP total

    lca_overall_table, lca_stream_table, lca_other_table, lca_unit_result, lca_final_table = lca_results()
    lca_unit_result_sums = np.abs(lca_unit_result).sum()
    lca_contributions = lca_unit_result.apply(lambda col: col*100 / lca_unit_result_sums[col.name], axis=0) # rescale the values so that the absolute value of the indicator sums to 100%
    lca_contributions = lca_contributions.sort_index(axis=1) # force the category names to be in the same order every time
    lca_contributions.columns = [
    'Photochemical Ox. Ecosystems',
    'Eutroph. Freshwater',
    'Ecotoxicity Freshwater',
    'Energy Resources',
    'Climate Change', 
    'Photochemical Ox. Human Health',
    'Human Toxicity Carc.',
    'Human Toxicity N-carc.',
    'Ionising Radiation',
    'Land Use',
    'Eutroph. Marine',
    'Ecotoxicity Marine',
    'Ozone Depletion',
    'Particulate Matter',
    'Meterial Resources',
    'Acidification Terrestrial',
    'Ecotoxicity Terrestrial',
    'Water Use'
    ]
    lca_contributions = lca_contributions.sort_index(axis=1, ascending=False)

    index_mapping = {'U1': 'Leaching', 'F1': 'Leaching', 'F2': 'Leaching',
        'P1': 'Concentration', 'F3': 'Concentration', 'H1': 'Concentration', 
        'P2': 'Refining', 'F4': 'Refining', 'H2': 'Refining',
        'RS': 'Separation','S1': 'Separation',
        'WT':'Wastewater Treatment', 'P3': 'Wastewater Treatment',
        'M1': 'Gypsum Credit',
        'M4': 'REO Credit',
        'M0': 'PG Credit'}
    # Replace the index values based on the mapping
    lca_contributions.index = lca_contributions.index.to_series().replace(index_mapping)

    # Group by index and aggregate 'Value' using sum
    indicator_contributions = lca_contributions.groupby(lca_contributions.index).agg('sum')

    # Add MSP to the other indicators
    testing_table = MSP_unit_table.copy() # make a copy of the MSP df (index: process sections, col: material cost, utility cost, gypsum credit, FOC, and capital recovery)
    testing_table.drop('Gypsum Credit', inplace=True, axis=1) # drop the gypsum credit column
    testing_table = np.abs(testing_table).sum(axis=1) # sum the across all the columns (reduce to one column)
    testing_table.loc['Gypsum Credit'] = -sys.get_market_value(fs_stream.gypsum)/product_mflow # add gypsum credit back as a row
    testing_table = testing_table.sort_index().multiply(100/np.abs(testing_table).sum()) # rescale the values so that the absolute value of the indicator sums to 100%
    indicator_contributions['Minimum Selling Price'] = testing_table # add the scaled MSP values to the scaled LCA values for plotting
    if fununit == 'PG':
        indicator_contributions = indicator_contributions.reindex(['Leaching', 'Concentration', 'Separation', 'Refining', 'Wastewater Treatment','Gypsum Credit', 'REO Credit', 'PG Credit']) # put the df in order so that it plots nicely
    elif fununit == 'Ln':
        indicator_contributions = indicator_contributions.reindex(['Leaching', 'Concentration', 'Separation', 'Refining', 'Wastewater Treatment','Gypsum Credit', 'PG Credit']) # put the df in order so that it plots nicely

    # Create a stacked bar chart
    # red f1777f
    # blue 60c1cf
    # green 79bf82
    # orange f98f60
    # purple a280b9
    # gray 90918e
    # yellow f3c354
    # black 403a48
    if fununit == 'PG':
        custom_colors = ['#f1777f', '#60c1cf', '#79bf82', '#f98f60', '#a280b9', '#90918e', '#966b6b']
    elif fununit == 'Ln':
        custom_colors = ['#f1777f', '#60c1cf', '#79bf82', '#f98f60', '#a280b9', '#90918e', '#403a48', '#98876e']

    aspect_ratio_LtoW = 1.25 # 6/10
    # Manually change the text font and size
    font = {'family': 'Calibri', 'size': 8}
    plt.rc('font', **font)

    fig, ax = plt.subplots(figsize=(width_one_col*1.25*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in)) # Matplotlib wants input in inches (width, length/height)
    indicator_contributions.T.plot(kind='barh', stacked=True, color=custom_colors, ax=ax)

    # Set labels and title
    ax.set_xlabel('Contribution to Indicator (%)')
    ax.set_xlim(-100,100)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=55, ha='right')
    ax.axvline(x=0, color='black', linewidth=0.8)
    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()
    #specify order of items in legend
    order = [0, 1, 2, 3, 4, 5, 6] # Controls the order of process sections in the figure legend. Each number corresponds to a process section
    #add legend to plot
    ax.legend().remove()  # This removes the default legend

    fig.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc=9, ncol = 2) # , bbox_to_anchor=(1.45, 0.75)
    fig.tight_layout() # rect=[0, 0, 0.95, 1]
    fig.subplots_adjust(top=0.8)
    # Show the plot
    fig.savefig(os.path.join(figures_path, f'Indicator_Contributions_{fununit}.tiff'), dpi=600)
    return MSP_table, MSP_unit_table, indicator_contributions