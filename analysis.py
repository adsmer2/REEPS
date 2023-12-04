"""
Rare Earth Element recovery from Phosphogypsum System (REEPS)

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""
# %%
import collections.abc
from collections.abc import Iterable
#hyper needs the four following aliases to be done manually.
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
collections.MutableSet = collections.abc.MutableSet
collections.MutableMapping = collections.abc.MutableMapping

# import packages
import qsdsan as qs
import biosteam as bst
from qsdsan.utils import (auom, clear_lca_registries)
from warnings import warn
import openpyxl
import os
import numpy as np 
import scipy as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# import functions from other files
from _component import *
from _lca_data import *
from model import *
from systems import *

from model_leaching import *
from bst_TEA_modified import * 

bw_path = os.path.dirname(__file__)
data_path = os.path.join(bw_path, 'data')

c_path = os.path.dirname(__file__)
figures_path = os.path.join(c_path, 'figures')

r_path = os.path.dirname(__file__)
results_path = os.path.join(r_path, 'results')

os.environ["PATH"] += os.pathsep + r'C:\Users\ajs8911\Miniconda3\envs\bw2\Library\bin' # (work) required to get graphviz system diagrams to save 
# os.environ["PATH"] += os.pathsep + r'C:\Users\adsme\miniconda3\envs\bw2\Library\bin' # (home) required to get graphviz system diagrams to save 

# Troubleshooting
# ------------------
# fs_stream.Ln_stream.show(flow = 'kg/hr') # shows the flow rates of components in the desired units
# fs_stream.Ln2O3.F_mass or F_mol or F_vol# shows total flow rate to many digits in kg/hr or kmol/hr or m3/hr

# %%
# =============================================================================
# Run the system and define what analysis to run
# =============================================================================

# fununit = 'PG'
# feedPG = 1000000
# REEcontent = 0.5/100
# num_ind_REEs = 9
# report = 'no' # "yes" or "no". Do you want to print results to excel? (saves 1 excel file to 'results' folder)
# num_samples = 1000
# uncertainty = 'no' # "yes" or "no". Do you want to make an kde-box/whisker plot? (saves 1 plot to 'figures' folder)
# sensitivity = 'no' # "yes" or "no". Do you want to make a bubble senstivity plot and parameter trend plots? (saves several plots to 'figures' folder)
# parameter = 'technological' # "contextual", "technological" or "all". sensitivity is split between these two types of parameters 
# optimization = 'no' # "yes" or "no". Do you want to make contour plots to pre-optimize the leaching unit? (saves 6 contour plots to 'figures' folder)
# desire_target ='no' # "yes" or "no". Do you want to identify if the technologies meet the targets? (saves 1 target plot to 'figures' folder)
# sys, lca, tea= create_system(fununit=fununit, feedPG=feedPG, REEcontent=REEcontent, num_ind_REEs=num_ind_REEs)
# flowsheet = qs.Flowsheet.flowsheet.default
# fs_stream = flowsheet.stream
# fs_unit = flowsheet.unit

def run_analysis(fununit, feedPG, REEcontent, num_ind_REEs,
                 report, num_samples, uncertainty, sensitivity, parameter, optimization, desire_target):
    fununit = fununit
    feedPG = feedPG
    REEcontent = REEcontent
    num_ind_REEs = num_ind_REEs
    sys, lca, tea= create_system(fununit=fununit, feedPG=feedPG, REEcontent=REEcontent, num_ind_REEs=num_ind_REEs)

    flowsheet = qs.Flowsheet.flowsheet.default
    fs_stream = flowsheet.stream
    fs_unit = flowsheet.unit

    # What are you looking to run?
    # Being able to choose yes or no for these helps reduce computation load for things you aren't looking for
    report = report # "yes" or "no". Do you want to print results to excel? (saves 1 excel file to 'results' folder)
    num_samples = num_samples
    uncertainty = uncertainty # "yes" or "no". Do you want to make an kde-box/whisker plot? (saves 1 plot to 'figures' folder)
    sensitivity = sensitivity # "yes" or "no". Do you want to make a bubble senstivity plot and parameter trend plots? (saves several plots to 'figures' folder)
    parameter = parameter # "contextual", "technological" or "all". sensitivity is split between these two types of parameters 
    optimization = optimization # "yes" or "no". Do you want to make contour plots to pre-optimize the leaching unit? (saves 6 contour plots to 'figures' folder)
    desire_target = desire_target # "yes" or "no". Do you want to identify if the technologies meet the targets? (saves 1 target plot to 'figures' folder)

    # =============================================================================
    # Print the system results to Excel
    # =============================================================================
    DataFrame = pd.DataFrame
    ExcelWriter = pd.ExcelWriter
    def tables_to_excel(tables, writer, sheet='Sheet1', n_row=1, row_spacing=2): 
        """
        Save a list of tables as an excel file and return the row number at which
        another consecutive table would start.
        
        Parameters
        ----------
        tables : iterable[pandas.DataFrame]
            Tables to save to excel.
        writer : pandas.ExcelWritter
            Writer to manage data stream to excel.
        sheet : str
            Name of sheet to save data.
        n_row : int
            Row number to begin saving data.
        row_spacing : int
            Number of rows between tables.
        
        Returns
        -------
        n_row : int
            Row number for next table.
            
        """
        row_spacing += 1 # Account for Python index offset
        for t in tables:
            label = t.columns.name
            t.to_excel(writer, sheet, 
                    startrow=n_row, index_label=label)
            n_row += len(t.index) + row_spacing
        return n_row

    def create_report(system, file='report.xlsx', dpi='300', tea=None, **stream_properties):
        """
        Save a system report as an xlsx file.
        
        Parameters
        ----------
        file : str
            File name to save report
        dpi : str, optional
            Resolution of the flowsheet. Defaults to '300'
        tea : TEA, optional
            Object for techno-economic analysis and cashflows. Defaults to the
            TEA object linked to the system.
        **stream_properties : str
            Additional stream properties and units as key-value pairs (e.g. T='degC', flow='gpm', H='kW', etc..)
            
        """
        writer = ExcelWriter(file)
        units = sorted(system.units, key=lambda x: x.line)
        cost_units = [i for i in units if i._design or i._cost]
        try:
            with bst.preferences.temporary() as p:
                p.reset()
                p.light_mode()
                system.diagram('thorough', file='flowsheet', dpi=str(dpi), format='png')
        except:
            diagram_completed = False
            warn(RuntimeWarning('failed to generate diagram through graphviz'), stacklevel=2)
        else:
            import PIL.Image
            try:
                # Assume openpyxl is used
                worksheet = writer.book.create_sheet('Flowsheet')
                flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
                worksheet.add_image(flowsheet, anchor='A1')
            except PIL.Image.DecompressionBombError:
                PIL.Image.MAX_IMAGE_PIXELS = int(1e9)
                flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
                worksheet.add_image(flowsheet, anchor='A1')
            except:
                # Assume xlsx writer is used
                try:
                    worksheet = writer.book.add_worksheet('Flowsheet')
                except:
                    warn("problem in saving flowsheet; please submit issue to BioSTEAM with"
                        "your current version of openpyxl and xlsx writer", RuntimeWarning)
                worksheet.insert_image('A1', 'flowsheet.png')
            diagram_completed = True
        
        def cost_table(tea): 
            """
            Return a cost table as a pandas DataFrame object.
            Parameters
            ----------
            units : iterable[Unit]
            Returns
            -------
            table : DataFrame
            """
            columns = ('Unit operation',
                    'Purchase cost (10^6 USD)',
                    'Utility Cost (10^6 USD/yr)',
                    'Additional OpEx Cost (10^6 USD/yr)',
                    'Total Operating Cost (10^6 USD/yr)')
            units = sorted([i for i in tea.system.units if i._design or i._cost], key=lambda x: x.line)
            operating_days = tea.operating_days
            N_units = len(units)
            array = np.empty((N_units, 5), dtype=object)
            IDs = []
            types = array[0:, 0]
            C_cap = array[0:, 1]
            C_ut = array[0:, 2]
            C_op = array[0:, 3]
            C_top = array[0:, 4]
            
            # Get data
            for i in range(N_units):
                unit = units[i]
                types[i] = unit.line
                C_cap[i] = unit.purchase_cost / 1e6
                C_ut[i] = unit.utility_cost * operating_days * 24  / 1e6
                C_op[i] = sum([i for i in unit.add_OPEX.values()]) * operating_days * 24  / 1e6
                C_top[i] = C_ut[i] + C_op[i]
                IDs.append(unit.ID)
            df = DataFrame(array, columns=columns, index=IDs)    
            if not tea.lang_factor:
                df['Installed cost (10^6 USD)'] = [u.installed_cost / 1e6 for u in units]
            return df
        
        # -------------------------------------------------------------------------------------------------------------------
        # Print to Excel
        # ---------------------------------------------------------------------------------------------------------------------

        # Stream tables
        # -----------------------------------------------------
        # Organize streams by chemicals first
        def _stream_key(s): # pragma: no coverage
            num = s.ID[1:]
            if num.isnumeric(): return int(num)
            else: return -1

        # Create a DataFrame with the total flowrate, source, sink, phase, and T of each stream
        def stream_table(streams, flow='kg/hr', percent=True, chemicals=None, **props):
            """
            Return a stream table as a pandas DataFrame object.

            Parameters
            ----------
            streams : array_like[Stream]
            flow : str
                Units for flow rate.
            props : str
                Additional stream properties and units as key-value pairs
            
            """
            
            # Prepare rows and columns
            ss = sorted(sorted([i for i in streams if i.ID], key=lambda i: i.ID), key=_stream_key)
            if not chemicals: 
                all_chemicals = tuple(set([i.chemicals for i in ss]))
                sizes = [(i, chemical.size) for i, chemical in enumerate(all_chemicals)]
                index, size = max(sizes, key=lambda x: x[1])
                chemicals = all_chemicals[index]
            n = len(ss)
            m = chemicals.size
            p = len(props)
            array = np.empty((m+p+5, n), dtype=object)
            IDs = n*[None]
            sources = array[0, :]
            sinks = array[1, :]
            phases = array[2, :]
            prop_molar_data = array[3:3+p+1,:]
            flows = array[p+3, :]
            array[p+4, :] = ''
            fracs = array[p+5:m+p+5, :]
            for j in range(n):
                s = ss[j]
                sources[j] = s.source.ID if s.source else '-'
                sinks[j] = s.sink.ID if s.sink else '-'
                IDs[j] = s.ID
                phase = ''
                for i in s.phase:
                    if i == 'l':
                        phase += 'liquid|'
                    elif i == 'L':
                        phase += 'LIQUID|'
                    elif i == 'g':
                        phase += 'gas|'
                    elif i == 's':
                        phase += 'solid|'
                phase = phase.rstrip('|')
                phases[j] = phase
                flow_j = s.get_flow(units=flow)
                flows[j] = net_j = flow_j.sum()
                if percent: net_j /= 100.
                fracs_j = flow_j/net_j if net_j > 1e-24 else 0
                if s.chemicals is chemicals:
                    fracs[:, j] = fracs_j
                else:
                    fracs[:, j] = 0.
                    fracs[chemicals.get_index(s.chemicals.IDs), j] = fracs_j
                i = 0
                for attr, units in props.items():
                    prop_molar_data[i, j] = s.get_property(attr, units)
                    i += 1
            index = (
                'Source', 
                'Sink',
                'Phase', 
                *[f'{attr} ({units})' for attr, units in props.items()], 
                f'flow ({flow})',
                ('Composition [%]:' if percent else 'Composition:'),
                *chemicals.IDs)
            return DataFrame(array, columns=IDs, index=index)

        # Now add the compositions of each of the streams [weight %]
        streams_by_chemicals = {}
        for i in system.streams:
            if not i: continue
            chemicals = i.chemicals
            if chemicals in streams_by_chemicals:
                streams_by_chemicals[chemicals].append(i)
            else:
                streams_by_chemicals[chemicals] = [i]
        stream_tables = []
        for chemicals, streams in streams_by_chemicals.items():
            stream_tables.append(stream_table(streams, chemicals=chemicals, T='K', **stream_properties))
        
        # Convert to mass flow rates as opposed to mass fractions
        columns = [e for sl in stream_tables for e in sl]

        index = []
        for j in np.array([i.chemicals for i in sys.streams][0]):
            index.append(j.ID)
        index

        mass_flows = np.reshape(np.array(stream_tables)[0,4,:], (len(np.array(stream_tables)[0,4,:]),1))

        result = []
        for row in np.arange(0,len(np.array(stream_tables)[0,:,0])):
            result.append(np.array(stream_tables)[0,row,:])
        result = pd.DataFrame(result)
        results = result.iloc[6:,:]

        mass_flow_result = pd.DataFrame(np.array(results)*mass_flows.T/100) # kg/hr. mass flowrates based on the above mass compositions and flow rates

        mass_flow_result.index = index 
        mass_flow_result.columns = columns

        # Complete the mass balance per each unit
        result = []
        result2 = []
        index = []
        for unit in fs_unit:
            result.append((np.sum([i.F_mass for i in unit.outs]) - np.sum([i.F_mass for i in unit.ins]))/np.sum([i.F_mass for i in unit.outs])*100)
            result2.append(np.sum([i.F_mass for i in unit.outs]) - np.sum([i.F_mass for i in unit.ins]))
            index.append(unit.ID)
        for i in np.arange(len(result)):
            if np.abs(result[i]) < 1E-3:
                result[i] = 0
            if np.abs(result2[i]) < 1E-6:
                result2[i] = 0
            else:
                pass
        mass_balance = pd.DataFrame(list(zip(index, result2, result))) # % difference in wt. Gives Output - Input/Output for each unit
        mass_balance.columns = ['Unit ID', 'difference (kg/hr)', 'wt pcnt difference (Out-In)']

        # Compile the tables for printing to Excel
        tables = [pd.DataFrame(stream_tables[0]), mass_flow_result, mass_balance]
        tables_to_excel(tables, writer, 'Stream Table')
        
        # General design requirements
        # ----------------------------------------------------------------
        def unit_result_tables(units): 
            """
            Return a list of results tables for each unit type.

            Parameters
            ----------
            units : iterable[Unit]
                
            Returns
            -------
            tables : list[DataFrame]
            
            """
            units = sorted(units, key=(lambda u: u.line))
            
            # Organize units by units of measure:
            organized = {}
            for u in units:
                uom = (*u._units.keys(), u.line)
                if uom in organized: organized[uom].append(u)
                else: organized[uom] = [u]
            
            # Make a list of tables, keeping all results with same keys in one table
            tables = []
            key = lambda u: u.ID
            for all_units in organized.values():
                # First table with units of measures
                all_units = sorted(all_units, key=key)
                u, *units = all_units
                key_hook = None
                table = u.results()
                table.columns.name = (u.line, '')
                tables.append(table)
            return tables

        results = unit_result_tables(cost_units)
        tables_to_excel(results, writer, 'Design requirements')

        # LCA Results
        # -------------------------------
        lca_stream_table, lca_other_table, lca_unit_result, lca_final_table = lca_results()
        tables = [lca_stream_table, lca_other_table, lca_unit_result, lca_final_table]
        tables_to_excel(tables, writer, 'LCA Results')

        # TEA Results (low level)
        # -------------------------------
        if tea is None: tea = system.TEA
        if tea:
            tea = system.TEA
            cost = cost_table(tea)
            stream_cost = stream_price_results()
            lst = [cost, stream_cost]
            # cost.to_excel(writer, 'Itemized costs')
            tables_to_excel(lst, writer, 'TEA low level')
            tea.get_cashflow_table().to_excel(writer, 'Cash flow')
        else:
            warn(f'Cannot find TEA object in {repr(system)}. Ignoring TEA sheets.',
                RuntimeWarning, stacklevel=2)

        # TEA Results (high level)
        # -------------------------------
        MSP_table, MSP_unit_table, indicator_contributions = MSP_results()
        tea_result = tea_results()
        lst = [tea_result, MSP_table, MSP_unit_table, indicator_contributions]
        tables_to_excel(lst, writer, 'Results high level')

        writer._save()
        if diagram_completed: os.remove("flowsheet.png")

    # -----------------------------------------------------------------------
    # LCA Results
    # -----------------------------------------------------------------------
    def power_utility_table(units): 
        """
        Return a pandas DataFrame object of power utilities.
        
        Parameters
        ----------
        units : iterable[Unit]
            
        """
        # Sort power utilities by unit type
        units = sorted(units, key=(lambda u: type(u).__name__))
        units = [u for u in units if u.power_utility]
        power_utilities = [u.power_utility for u in units]
        length = len(power_utilities)
        data = []
        for i, u, pu in zip(range(length), units, power_utilities):
            data.append((u.line, pu.rate, pu.cost))
        return DataFrame(data, index=[u.ID for u in units if u.power_utility],
                        columns=('Unit Operation', 'Rate (kW)', 'Cost (USD/hr)'))

    def lca_results():
        # get impacts for each stream and utility
        lca_stream_table = lca.get_impact_table('Stream') # get lca indicator results for each flow in/out of the system
        lca_other_table = lca.get_impact_table('Other') # get lca indicator results for utilities being used in the system
        
        # Get the impacts for each unit in the system
        results = []
        for u in fs_unit: # for each unit in the system, get the lca results and save it to results
            result = lca.get_unit_impacts(units=u)
            myKeys = list(result.keys())
            myKeys.sort()
            sorted_dict = {i: result[i] for i in myKeys}
            myValues = sorted_dict.values()
            d1 = pd.DataFrame(myValues, index = myKeys)
            result = d1.values.tolist()
            results.append(result)
        lca_unit_result = pd.DataFrame(np.reshape(np.array(results),(len([u for u in fs_unit]), len(lca.indicators)))) # Make results a dataframe
        # These lca results from a preexisting function do not consider the utilities unit specifically and instead incorrectly add the electricity and heating impacts to each unit (needs to be further processed)
        col_names = [i.category for i in lca.indicators] # Get the indicator names for the columns of the df
        col_names = [category.title() for category in col_names] # Capitalize the indicator names
        col_names.sort()
        ind_names = [i.ID for i in fs_unit] # get unit ID names for the rows of the df
        lca_unit_result.columns = col_names
        lca_unit_result.index = ind_names

        # get the electricity consumption for each unit in the system
        df_util = power_utility_table(fs_unit).iloc[:,1]
        elec_unit_lca = pd.DataFrame(df_util/np.sum(df_util)) # dataframe with the fraction of the total electricity being used by each unit

        # get the natural gas consumption for each unit
        heatNG_index = ['U1', 'H1', 'H2']
        heatNG_sum = np.sum([fs_unit.U1.heat_duty, fs_unit.H1.heat_duty, fs_unit.H2.heat_duty])
        heatNG_values = [fs_unit.U1.heat_duty/heatNG_sum, fs_unit.H1.heat_duty/heatNG_sum, fs_unit.H2.heat_duty/heatNG_sum] # get the fraction of the total natural gas being used
        heatNG_unit_lca = pd.DataFrame(heatNG_values) # dataframe with the fraction of the total natural gas being used by each unit
        heatNG_unit_lca.index = heatNG_index

        # Get the total value of impact from electricity and natural gas
        df = lca.get_impact_table('Other')
        df2 = df.loc[:,['[' in i for i in df.columns]] # get only the total value of impact (not the ratio of impact from electricity and NG) for each impact category
        dfe = df2.loc[['electricity' in i for i in df.index], :] # total impact of electricity consumption
        dfNG = df2.loc[['NG' in i for i in df.index], :] # total impact of natural gas heating
        dfSum = df2.loc[['Sum' in i for i in df.index], :] # total impact of utilities (heating + electricity)

        # multiply the total impact of electricity times the fraction of electricity used for each unit
        df_elec = pd.DataFrame(np.array(elec_unit_lca).reshape(len(elec_unit_lca),1)*np.array(dfe)) # df of the impact from electricity consumption of each unit
        df_elec.index = elec_unit_lca.index
        df_elec.columns = dfe.columns
        df_elec = df_elec.sort_index(axis=1)

        # multiply the total impact of natural gas for heat times the fraction of NG used for each unit
        df_NG = pd.DataFrame(np.array(heatNG_unit_lca)*np.array(dfNG)) # df of the impact from natural gas heating of each unit
        df_NG.index = heatNG_unit_lca.index
        df_NG.columns = dfNG.columns
        df_NG = df_NG.sort_index(axis=1)

        # subtract the sum of impacts from utilities from the impact of each unit (correct the existing function's incorrect addition of utility impacts)
        lca_unit_result = lca_unit_result - dfSum.sort_index(axis=1).values 

        # Now add the utility impacts back correctly for each unit specifically
        for row in lca_unit_result.index: # for units that use electricity, add the impact of that electricity use to the system's overall impact for each impact category
            for row2 in df_elec.index:
                if row == row2:
                    lca_unit_result.loc[row,:] = lca_unit_result.loc[row,:].values.reshape(1,len(lca_unit_result.columns)) + df_elec.loc[row2,:].values.reshape(1,len(df_elec.columns))
        for row in lca_unit_result.index: # for units that use natural gas for heat, add the impact of that natural gas use to the system's overall impact for each impact category
            for row3 in df_NG.index:
                if row == row3:
                    lca_unit_result.loc[row,:] = lca_unit_result.loc[row,:].values.reshape(1,len(lca_unit_result.columns)) + df_NG.loc[row3,:].values.reshape(1,len(df_NG.columns))
        # ***double checked this fix by adding the sum of stream impacts plus the sum of utility impacts and comparing to the sum of the unit impact results for each category***

        # make a final summary table for lca results
        # -------------------------------------------
        lca_ind = ['PG Feed Flow', 'Ln2O3 Product Flow','Total GW EI99/PG (compare kulczycka 1.56 pt)', 'Total GW ReCiPe/PG', 'Total GW ReCiPe/Ln']
        lca_value = [fs_stream.rawPG.F_mass, fs_stream.Ln2O3.F_mass, lca.total_impacts['A_ClimateChange']/(fs_stream.rawPG.F_mass*tea.operating_days*24)*1000, lca.total_impacts['GWP100']/(fs_stream.rawPG.F_mass*tea.operating_days*24), lca.total_impacts['GWP100']/(fs_stream.Ln2O3.F_mass*tea.operating_days*24)]
        lca_unit = ['kg/hr', 'kg/hr', 'EI99 pts/mt PG. Note: not accurate unless the correct functional unit is chosen', 'kg CO2-eq/kg PG. Note: not accurate unless the correct functional unit is chosen', 'kg CO2-eq/kg Ln. Note: not accurate unless the correct functional unit is chosen']
        lca_final_table = pd.DataFrame(list(zip(lca_value, lca_unit)))
        lca_final_table.columns = ['Value','Unit']
        lca_final_table.index = lca_ind

        return lca_stream_table, lca_other_table, lca_unit_result, lca_final_table

    # -----------------------------------------------------------------------
    # TEA Results (low level)
    # -----------------------------------------------------------------------
    def stream_price_results():
        ind_names = [i.linked_stream.ID for i in lca.stream_inventory]
        prices = [i.price for i in lca.stream_inventory]
        streamFlows = [i.linked_stream.F_mass for i in lca.stream_inventory]
        streamCost = np.array(prices)*np.array(streamFlows)*(tea.duration[1]-tea.duration[0])*tea.operating_hours
        streamCost_annual = streamCost/(tea.duration[1]-tea.duration[0])
        pResult = prices + streamFlows + list(streamCost)
        stream_cost_result = pd.DataFrame(list(zip(prices, streamFlows, streamCost, streamCost_annual)))
        col_names = ['Price $/kg', 'Mass Flow (kg/hr)', 'Flow Value (total $ over duration of TEA)', 'Flow Value ($/year)']
        ind_names = [i.linked_stream.ID for i in lca.stream_inventory]
        stream_cost_result.columns = col_names
        stream_cost_result.index = ind_names
        stream_cost_result
        
        return stream_cost_result

    # -----------------------------------------------------------------------
    # TEA Results (High level)
    # -----------------------------------------------------------------------
    def tea_results():
        tea_ind = ['TCI','Sales','FOC','VOC','AOC','Net Earnings','ROI','PBP','IRR','DCFRR','NPV','MSP']
        tea_value = [tea.TCI/1e6, tea.sales/1e6, tea.FOC/1e6, tea.VOC/1e6, tea.AOC/1e6, tea.net_earnings/1e6, tea.ROI, tea.PBP, tea.IRR, tea.solve_IRR(), tea.NPV/1e6, tea.solve_price(fs_stream.Ln2O3)]
        tea_unit = ['million $','million $/year','million $/year','million $/year','million $/year','million $/year','1/year','years','-','-','million $','$/kg']
        tea_table = pd.DataFrame(list(zip(tea_value, tea_unit)))
        tea_table.columns = ['Value','Unit']
        tea_table.index = tea_ind
        return tea_table

    def MSP_results():
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
        aspect_ratio_LtoW = 6/4
        cm_to_in = 1/2.54  # centimeters in inches
        width_one_col = 8.3 # cm. Width for a one column figure
        width_two_col = 17.1 # cm. Width for a two column figure
        max_length = 23.3 # cm. The maximum lenght a figure can be

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))
        ax = MSP_table.drop('MSP', axis=0).T.plot(kind='bar', stacked=True, color=custom_colors, ax=ax)
        # Set labels and title
        ax.set_ylabel('Minimum Selling Price (USD/kg REO)')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1) ) # , ncol=len(combined_df.index)
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
        # gray 403a48
        custom_colors = ['#82cfd0', '#98876e', '#403a48', '#fcb813', '#3ba459']

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

        lca_stream_table, lca_other_table, lca_unit_result, lca_final_table = lca_results()
        lca_unit_result_sums = np.abs(lca_unit_result).sum()
        lca_unit_result.drop('Human Health', inplace=True, axis=1) # Gets rid of duplicate EI-99 LCIA result for GW
        lca_contributions = lca_unit_result.apply(lambda col: col*100 / lca_unit_result_sums[col.name], axis=0) # rescale the values so that the absolute value of the indicator sums to 100%

        index_mapping = {'U1': 'Leaching', 'F1': 'Leaching', 'F2': 'Leaching',
            'P1': 'Concentration', 'F3': 'Concentration', 'H1': 'Concentration', 
            'P2': 'Refining', 'F4': 'Refining', 'H2': 'Refining',
            'RS': 'Separation','S1': 'Separation',
            'WT':'Wastewater Treatment', 'P3': 'Wastewater Treatment',
            'M1': 'Gypsum Credit',
            'M4': 'REO Credit',
            'M0': 'PG Remediation Credit'}
        # Replace the index values based on the mapping
        lca_contributions.index = lca_contributions.index.to_series().replace(index_mapping)
        lca_contributions['Natural Land Transformation'] = lca_contributions.loc[:, 'Natural Land Transformation'].multiply(-1) # for some reason CFs for NLTP are inverted (*-1). This appears to be an ecoinvent 3.8 LCIA method list issue not a brightway2 issue.

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
            indicator_contributions = indicator_contributions.reindex(['Leaching', 'Concentration', 'Separation', 'Refining', 'Wastewater Treatment','Gypsum Credit', 'REO Credit']) # put the df in order so that it plots nicely
        else:
            indicator_contributions = indicator_contributions.reindex(['Leaching', 'Concentration', 'Separation', 'Refining', 'Wastewater Treatment','Gypsum Credit', 'PG Remediation Credit']) # put the df in order so that it plots nicely

        # Create a stacked bar chart
        # red f1777f
        # blue 60c1cf
        # green 79bf82
        # orange f98f60
        # purple a280b9
        # gray 90918e
        # yellow f3c354
        # black 403a48
        custom_colors = ['#f1777f', '#60c1cf', '#79bf82', '#f98f60', '#a280b9', '#90918e', '#403a48', '#403a48']

        aspect_ratio_LtoW = 1 # 6/10
        fig, ax = plt.subplots(figsize=(width_two_col*cm_to_in, width_two_col*aspect_ratio_LtoW*cm_to_in)) # Matplotlib wants input in inches (width, length/height)
        ax = indicator_contributions.T.plot(kind='bar', stacked=True, color=custom_colors, ax=ax)

        # Set labels and title
        ax.set_ylabel('Contribution to Indicator (%)')
        ax.set_ylim(-100,100)
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=55, ha='right')
        ax.axhline(y=0, color='black', linewidth=0.8)
        #get handles and labels
        handles, labels = plt.gca().get_legend_handles_labels()
        #specify order of items in legend
        order = [0, 1, 2, 3, 4, 5, 6] # Controls the order of process sections in the figure legend. Each number corresponds to a process section
        #add legend to plot
        ax.legend().remove()  # This removes the default legend
        fig.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc=9, ncol = 4) # , bbox_to_anchor=(1.45, 0.75)
        fig.tight_layout() # rect=[0, 0, 0.95, 1]
        fig.subplots_adjust(top=0.9)
        # Show the plot
        fig.savefig(os.path.join(figures_path, f'Indicator_Contributions_{fununit}.tiff'), dpi=600)
        return MSP_table, MSP_unit_table, indicator_contributions

    if report == 'yes':
        create_report(system=sys, file=os.path.join(results_path, f'{sys.ID}_results_{fununit}_{feedPG}_{REEcontent}_{num_ind_REEs}.xlsx'))


    # ===========================================================================================================
    # Uncertainty and Sensitivity for base model
    # ===========================================================================================================

    target = 'no' # To let the first model run without considering the additional uncertainty of considering a larger range of the target parameter
    def run_uncertainty(sys, fununit):
        # Create model
        # -------------------------
        model_uncertainty = create_model(sys, fununit, 'all', target) # set system, functional unit, and type of parameter for analysis

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
        ydata = model_uncertainty.table.loc[:,('LCA', 'Global Warming [kg CO2-Eq/kg PG]')]
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

    def run_sensitivity(sys, fununit):
        # Sensitivity Analysis
        # -------------------------
        # use above model to get correlations 
        model_sensitivity = create_model(sys, fununit, parameter, 'no') # set system, functional unit, and type of parameter for analysis

        np.random.seed(3221) # setting the seed ensures you get the same sample

        samples = model_sensitivity.sample(N=num_samples, rule='L')
        model_sensitivity.load_samples(samples)
        model_sensitivity.evaluate()
        r_df, p_df = qs.stats.get_correlations(model_sensitivity, kind='Spearman')

        # # sort the parameters by alphabetically by unit ID then parameter name
        # lst1 = [i.element for i in model_sensitivity.get_parameters()]
        # lst2 = [i.name for i in model_sensitivity.get_parameters()]
        # key_parameters = sorted(list(zip(lst1, lst2)), key=lambda x: x[0])

        # Filter out parameters that only meet a certain threshold
        def filter_parameters(model, df, threshold):
            new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
            filtered = new_df.dropna(how='all')
            param_dct = {p.name_with_units:p for p in model.get_parameters()}
            parameters = set(param_dct[i[1]] for i in filtered.index)
            return list(parameters)
        key_parameters2 = filter_parameters(model_sensitivity, r_df, threshold=0.3) # Only want parameters with Spearman's rho >= 0.4 or <= -0.4

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

        filtered_unit_name = [i.element for i in key_parameters2]
        filtered_param_name = [i.name for i in key_parameters2]
        key_parameters2 = sorted(list(zip(filtered_unit_name, filtered_param_name)), key=lambda x: x[0])
        param_names = _update_input(np.array(key_parameters2)[:,1], df.index)
        # param_names = [i.name for i in key_parameters2]
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

        # make DataFrame
        data = {'Sign': corr_df['Sign']}
        df_gpt = pd.DataFrame(data)

        # Function to categorize values
        def categorize_sign(value):
            if value < 0:
                return '$-$'
            elif value > 0:
                return '$+$'
            else:
                return 'Zero'

        # Apply the function to create a new column
        corr_df['Sign'] = df_gpt['Sign'].apply(categorize_sign)

        # Remove all metrics that aren't economic because enviornmental impacts are unchanged by price changes
        if parameter == 'contextual':
            corr_df = corr_df[~corr_df['metric'].isin(metric_names[3:])] 

        # Begin plotting
        def _plot_corr_bubble(corr_df, ratio, **kwargs):
            plt.style.use('default')

            if parameter == 'technological':
                margin_x = kwargs['margin_x'] if 'margin_x' in kwargs.keys() else 0.05
                margin_y = kwargs['margin_y'] if 'margin_y' in kwargs.keys() else 0.05
                kwargs = {i: kwargs[i] for i in kwargs.keys() if 'margin' not in i}

                keys = ('height', 'aspect', 'sizes', 'size_norm', 'edgecolor') # , 'palette'
                values = (9+ratio, 1, (0, 1000), (0, 2.5), '0.5') # , ['#60c1cf', '#f1777f']
            elif parameter == 'contextual':
                margin_x = kwargs['margin_x'] if 'margin_x' in kwargs.keys() else 0.2/ratio
                margin_y = kwargs['margin_y'] if 'margin_y' in kwargs.keys() else 0.05
                kwargs = {i: kwargs[i] for i in kwargs.keys() if 'margin' not in i}

                keys = ('height', 'aspect', 'sizes', 'size_norm', 'edgecolor') # , 'palette'
                values = (7+ratio, 0.7, (0, 1000), (0, 2.5), '0.5') # , ['#60c1cf', '#f1777f']
            else: 
                RuntimeError(f'parameter={parameter} is not "technological" or "contextual". Please define as one of these two.')

            for num, k in enumerate(keys):
                kwargs.setdefault(keys[num], values[num])

            g = sns.relplot(data=corr_df, x='metric', y='parameter',
                            hue='Sign', size='Correlation', palette= {'$+$':'#60c1cf', '$-$':'#f1777f'}, **kwargs)

            g.set(xlabel='', ylabel='', aspect=1)
            g.ax.margins(x=margin_x, y=margin_y)

            for label in g.ax.get_xticklabels():
                label.set_rotation(90)

            for artist in g.legend.legendHandles:
                artist.set_edgecolor('1')

            for key in g.ax.spines.keys():
                g.ax.spines[key].set(color='k', linewidth=0.5, visible=True)

            g.ax.grid(True, which='major', color='k',linestyle='--', linewidth=0.3)
            g.tight_layout()
            
            if parameter == 'technological':
                sns.move_legend(g, 'center right', bbox_to_anchor=(1.025, 0.55))
            elif parameter == 'contextual':
                sns.move_legend(g, 'center right', bbox_to_anchor=(0.85, 0.55))
            else: 
                RuntimeError(f'parameter={parameter} is not "technological" or "contextual". Please define as one of these two.')

            return g

        # g = _plot_corr_bubble(corr_df, len(metric_names)/len(param_names))
        # g.savefig(os.path.join(figures_path, f'bubble_sensitivity_{fununit}_{parameter}.tiff'), dpi=600)

        pivot_corr_df = corr_df2.pivot(index="parameter", columns="metric", values="Correlation")
        pivot_corr_df = pivot_corr_df.sort_index(key=lambda x: [i.split('(')[1] for i in x])
        pivot_corr_df = pivot_corr_df.reindex(columns=['Agricultural Land Occupation',
        'Fossil Depletion',
        'Freshwater Ecotoxicity',
        'Freshwater Eutrophication',
        'Global Warming',
        'Human Toxicity',
        'Ionising Radiation',
        'Marine Ecotoxicity',
        'Marine Eutrophication',
        'Metal Depletion',
        'Natural Land Transformation',
        'Ozone Depletion',
        'Particulate Matter Formation',
        'Photochemical Oxidant Formation',
        'Terrestrial Acidification',
        'Terrestrial Ecotoxicity',
        'Urban Land Occupation',
        'Water Depletion',
        'NPV15',
        'IRR',
        'MSP'])

        # Create figures the correct size for publication
        aspect_ratio_LtoW = 0.65
        cm_to_in = 1/2.54  # centimeters in inches
        width_one_col = 8.3 # cm. Width for a one column figure
        width_two_col = 17.1 # cm. Width for a two column figure
        max_length = 23.3 # cm. The maximum lenght a figure can be

        fig, ax = plt.subplots(figsize=(width_two_col*cm_to_in, width_two_col*aspect_ratio_LtoW*cm_to_in))
        sns.heatmap(data=pivot_corr_df, cmap=sns.color_palette("vlag_r", as_cmap=True))
        ax.set(xlabel='', ylabel='')
        fig.tight_layout()
        fig.savefig(os.path.join(figures_path, f'Sensitivity_contributions_heatmap.tiff'), dpi=600)

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
        return model_sensitivity

    def ind_trend_analysis(model, p_name):
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
        # plt.show()


    if uncertainty == 'yes':
        model_uncertainty = run_uncertainty(sys, fununit)
        qs.Model._reset_system(model_uncertainty)
        qs.Model.metrics_at_baseline(model_uncertainty)
    if sensitivity == 'yes':
        model_sensitivity = run_sensitivity(sys, fununit)
        if parameter == 'technological':
            ind_trend_analysis(model_sensitivity, 'Acid Concentration (U1)') # input parameter name as it appears in model.py @param()
            ind_trend_analysis(model_sensitivity, 'Solvent to Solid Ratio (U1)')
            ind_trend_analysis(model_sensitivity, 'Sodium Hydroxide Feed (P3)')
        qs.Model._reset_system(model_sensitivity)
        qs.Model.metrics_at_baseline(model_sensitivity)


    # =============================================================================
    # Preoptimization Contour Plots
    # =============================================================================

    # def leaching_solventRatio_temp_NPV(sys, fununit):
    #     model = create_model_leaching(sys, fununit)
    #     # order of model parameter output [time, temp, conc, solventRatio]

    #     y = np.linspace(25, 70, 20) # y axis variable (temperature)
    #     x = np.linspace(2, 7, num=len(y)) # y axis variable (time)
    #     f1 = [5]*len(y) # fixed parameter not changing
    #     f2 = [120]*len(y) # fixed parameter not changing
    #     result = []
    #     i = 0

    #     while i < len(y):
    #         xi = [x[i] for n in range(len(y))]
    #         samples = pd.DataFrame(columns = ['x', 'y', 'f1', 'f2'])
    #         samples['x'], samples['y'], samples['f1'], samples['f2'] = f2, y, f1, xi # model takes columns in order initiated in model.py
    #         model.load_samples(samples.to_numpy())
    #         model.evaluate()
    #         result.extend(model.table.loc[:,('TEA', 'NPV [MM USD/kg]')].values.tolist())
    #         i += 1

    #     X, Y = np.meshgrid(x, y)
    #     Z = np.array(result).reshape(len(y),len(y)).T

    #     fig, ax = plt.subplots()
    #     cp = plt.contourf(X,Y,Z, 20, cmap='inferno')
    #     fig.colorbar(cp, label=f'NPV (MM USD)')
    #     ax.set(xlabel='Liquid/Solid Ratio', ylabel='Leaching Temperature (deg C)')
    #     fig.savefig(os.path.join(figures_path, f'test_leaching_temp_solventRatio_NPV_{fununit}.png'), dpi=300)

    # def leaching_solventRatio_time_NPV(sys, fununit):
    #     model = create_model_leaching(sys, fununit)
    #     # order of model parameter output [time, temp, conc, solventRatio]

    #     y = np.linspace(25, 205, 20) # y axis variable (temperature)
    #     x = np.linspace(2, 7, num=len(y)) # y axis variable (time)
    #     f1 = [5]*len(y) # fixed parameter not changing
    #     f2 = [50]*len(y) # fixed parameter not changing
    #     result = []
    #     i = 0

    #     while i < len(y):
    #         xi = [x[i] for n in range(len(y))]
    #         samples = pd.DataFrame(columns = ['x', 'y', 'f1', 'f2'])
    #         samples['x'], samples['y'], samples['f1'], samples['f2'] = y, f2, f1, xi # model takes columns in order initiated in model.py
    #         model.load_samples(samples.to_numpy())
    #         model.evaluate()
    #         result.extend(model.table.loc[:,('TEA', 'NPV [MM USD/kg]')].values.tolist())
    #         i += 1

    #     X, Y = np.meshgrid(x, y)
    #     Z = np.array(result).reshape(len(y),len(y)).T

    #     fig, ax = plt.subplots()
    #     cp = plt.contourf(X,Y,Z, 20, cmap='inferno')
    #     fig.colorbar(cp, label=f'NPV (MM USD)')
    #     ax.set(xlabel='Liquid/Solid Ratio', ylabel='Leaching Time (min)')
    #     fig.savefig(os.path.join(figures_path, f'test_leaching_solventRatio_time_NPV_{fununit}.png'), dpi=300)

    # # Old calculation
    # leaching_solventRatio_temp_NPV(sys, fununit)
    # leaching_solventRatio_time_NPV(sys, fununit)

    def optimization_leaching(system, fununit, xdata, ydata, f1data, f2data, indicator, nsamples): # set system and functional unit
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

        model = create_model_leaching(sys, fununit) # order of model parameter output [time, temp, conc, solventRatio]
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
        fig, ax = plt.subplots()
        cp = plt.contourf(X,Y,Z, 20, cmap=ind_cmap)
        fig.colorbar(cp, label=ind_axis_label)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        fig.savefig(os.path.join(figures_path, f'leaching_{xdata}_{ydata}_{indicator}_{fununit}.tiff'), dpi=600)

        qs.Model._reset_system(model)
        qs.Model.metrics_at_baseline(model)

    if optimization == 'yes':
        # time   temp   conc   solventRatio   NPV   GWP
        # Optimization by NPV
        optimization_leaching(system=sys, fununit=fununit, xdata='solventRatio', ydata='temp', f1data='conc', f2data='time', indicator='NPV', nsamples=20)
        optimization_leaching(system=sys, fununit=fununit, xdata='solventRatio', ydata='time', f1data='conc', f2data='temp', indicator='NPV', nsamples=20)
        optimization_leaching(system=sys, fununit=fununit, xdata='solventRatio', ydata='conc', f1data='time', f2data='temp', indicator='NPV', nsamples=20)
        optimization_leaching(system=sys, fununit=fununit, xdata='conc', ydata='time', f1data='temp', f2data='solventRatio', indicator='NPV', nsamples=20)
        optimization_leaching(system=sys, fununit=fununit, xdata='conc', ydata='temp', f1data='time', f2data='solventRatio', indicator='NPV', nsamples=20)
        optimization_leaching(system=sys, fununit=fununit, xdata='temp', ydata='time', f1data='conc', f2data='solventRatio', indicator='NPV', nsamples=20)


    # =============================================================================
    # Target Analysis
    # =============================================================================
    if desire_target == 'yes':
        
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


        # # Target 2
        # # --------------------------
        # sys, lca, tea= create_system(fununit=fununit, feedPG=feedPG, REEcontent=REEcontent, num_ind_REEs=num_ind_REEs)
        # target = 'REE Recovery (S1)'
        # model = create_model(sys, fununit, 'all', target=target) # set system, functional unit, and type of parameter for analysis

        # np.random.seed(3221) # setting the seed ensures you get the same sample

        # samples = model.sample(N=num_samples, rule='L')
        # model.load_samples(samples)
        # model.evaluate()

        # # Define parameter of interest
        # p_name = target
        # ind_parameter = [i.name for i in model._parameters].index(p_name)
        # p_units = model._parameters[ind_parameter].units
        
        # # Define range of parameter
        # lower = model._parameters[ind_parameter].distribution.lower[0]
        # upper = model._parameters[ind_parameter].distribution.upper[0]
        # num_bins = 65
        
        # # Define the bin width and make a list of bin centers
        # bin_width = (upper-lower)/num_bins
        # bin_centers = np.linspace(lower + bin_width/2, upper - bin_width/2, num=int((upper - lower) / bin_width))
        # mean_output_parameter_binned = []
        # std_output_parameter_binned = []

        # # Gather the mean and standard deviation for each bin from the model results
        # for bin_center in bin_centers:
        #     bin_start = bin_center - bin_width/2
        #     bin_end = bin_center + bin_width/2
        #     bin_indices = np.logical_and(model.table.loc[:, (model._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values >= bin_start, 
        #                                 model.table.loc[:, (model._parameters[ind_parameter].element, f'{p_name} [{p_units}]')].values < bin_end) # Gather indices of the parameter values within the bin
        #     bin_data = model.table.loc[:,('TEA', 'NPV15 [MM USD/kg]')].values[bin_indices] # using the relevant bin indices, create a list of indicator values.  .loc[:,('TEA', 'IRR [%]')]
        #     mean_output_parameter_binned.append(np.mean(bin_data)) # take the mean of the indicator values for this bin
        #     std_output_parameter_binned.append(np.std(bin_data)) # take the SD of the indicator values for this bin
        # # Smooth the data using cubic spline interpolation 
        # # Note: interp fills the gap between 'lower'/smallest 'bin_center' and 'upper'/higher 'bin_center' leading to a small 'flatline' at either end of plot. Can be minimized with more bins.
        # smooth_x_values = np.linspace(lower, upper, num=1000)  # Increase the number of points for smooth curves
        # mean_interpolator = np.interp(smooth_x_values, bin_centers, mean_output_parameter_binned)
        # std_interpolator = np.interp(smooth_x_values, bin_centers, std_output_parameter_binned)

        # # Create the line plot with uncertainty using plot and fill_between
        # plt.style.use('default')
        # fig, ax = plt.subplots()
        # ax.plot(smooth_x_values, mean_interpolator, color='#f98f60', label='Mean NPV15')
        # ax.fill_between(smooth_x_values, mean_interpolator - 2*std_interpolator, mean_interpolator + 2*std_interpolator,color='#f98f60', alpha=0.3, label='Confidence Interval (95%)')
        # # Add labels and title
        # ax.set(xlabel=f'{p_name} ({p_units})', ylabel='NPV15 (MM USD)', xlim=(lower,upper))
        # ax.grid(visible=True)
        # ax.ticklabel_format(axis='x', style='sci')
        # # Show the legend
        # fig.tight_layout()
        # ax.legend(loc='upper left')
        # # Save the plot
        # fig.savefig(os.path.join(figures_path, f'Target_{p_name}.tiff'), dpi=600)

    return feedPG, REEcontent, tea.solve_price(fs_stream.Ln2O3)


run_analysis(fununit='PG', feedPG=1000000, REEcontent=0.5/100, num_ind_REEs=9,
             report='yes', num_samples=500, uncertainty='no', sensitivity='no', parameter='technological', optimization='no', desire_target='no')

# # Scenario Analysis
# x = np.linspace(0.02/100, 0.9/100, 20) # wt fraction. REE content
# y = np.linspace(0.1*1e6, 2*1e6, 20) # M kg/hr. capacity

# resultCapacity = []
# resultREEContent = []
# resultMSP = []
# for i in x:
#     for j in y:
#         dataCapacity, dataREEContent, dataMSP = run_analysis(fununit='PG', feedPG= j, REEcontent= i, num_ind_REEs=9,
#              report='no', num_samples=1, uncertainty='no', sensitivity='no', parameter='technological', optimization='no', desire_target='no')
#         resultCapacity.append(dataCapacity)
#         resultREEContent.append(dataREEContent)
#         resultMSP.append(dataMSP)

# X, Y = np.meshgrid(x*100, y/1e6)
# Z = np.array(resultMSP).reshape(len(y),len(y)).T

# plt.style.use('default')
# aspect_ratio_LtoW = 2880/3840
# cm_to_in = 1/2.54  # centimeters in inches
# width_one_col = 8.3 # cm. Width for a one column figure
# width_two_col = 17.1 # cm. Width for a two column figure
# max_length = 23.3 # cm. The maximum lenght a figure can be
# fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))
# cp = plt.contourf(X,Y,Z, levels=np.linspace(0,400,30), cmap='inferno_r', extend='max')
# fig.colorbar(cp, label=f'MSP ($/kg REO)', ticks=np.arange(0,370,60)) 

# cp2 = plt.contour(X,Y,Z, levels=[51.5, 51.5*2], colors='black', linestyles='dashed')
# fmt = {}
# strs = [' MSP  51.5 ', ' MSP  103 '] 
# for l, s in zip(cp2.levels, strs):
#     fmt[l] = s
# ax.clabel(cp2, inline=True, fmt=fmt)

# ax.set(xlabel='REE Content of the PG (wt %)', ylabel='Capacity (M kg PG/hr)')
# fig.tight_layout()
# fig.savefig(os.path.join(figures_path, f'Scenario_analysis.tiff'), dpi=600)
