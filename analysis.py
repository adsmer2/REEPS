"""
Rare Earth Element recovery from Phosphogypsum System (REEPS)

The Pennsylvania State University
Chemical Engineering Department
S2D2 Lab (Dr. Rui Shi)
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
from matplotlib.ticker import ScalarFormatter

# import essential functions from other files
from _component import *
from _lca_data import *
from model import *
from systems import *
from model_leaching import *
from bst_TEA_modified import * 

# import the analysis functions from other files
from analysis_MSP_contributions import *
from analysis_uncertainty import *
from analysis_sensitivity import *
from analysis_indicator_trend import *
from analysis_optimization_leaching import *
from analysis_target import *
from analysis_scenario import *

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
# num_samples = 400
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
                 report, num_samples, uncertainty, sensitivity, parameter, optimization, desire_target, desire_scenario):
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
        # MSP_table, MSP_unit_table, indicator_contributions = MSP_results()
        MSP_table, MSP_unit_table, indicator_contributions = analysis_MSP_contributions(sys=sys, tea=tea, fs_stream=fs_stream, fs_unit=fs_unit, fununit=fununit, lca_results=lca_results, figures_path=figures_path)
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

    # Create analysis excel report and some baseline figures
    if report == 'yes':
        create_report(system=sys, file=os.path.join(results_path, f'{sys.ID}_results_{fununit}_{feedPG}_{REEcontent}_{num_ind_REEs}.xlsx'))


    # Uncertainty Analysis
    # -------------------------
    if uncertainty == 'yes':
        model_uncertainty = analysis_uncertainty(sys, fununit, num_samples, figures_path)
        qs.Model._reset_system(model_uncertainty)

    # Sensitivity Analysis
    # -------------------------
    if sensitivity == 'yes':
        model_sensitivity = analysis_sensitivity(sys, fununit, parameter, num_samples, figures_path)
        if parameter == 'technological':
            analysis_indicator_trend(model_sensitivity, 'Acid Concentration (U1)', figures_path=figures_path) # input parameter name as it appears in model.py @param()
            analysis_indicator_trend(model_sensitivity, 'Solvent to Solid Ratio (U1)', figures_path=figures_path)
            analysis_indicator_trend(model_sensitivity, 'Sodium Hydroxide Feed (P3)', figures_path=figures_path)
        qs.Model._reset_system(model_sensitivity)


    # Optimization Contour Plots
    # -----------------------------
    if optimization == 'yes':
        # time   temp   conc   solventRatio   NPV   GWP
        # Optimization by NPV
        analysis_optimization_leaching(system=sys, fununit=fununit, xdata='solventRatio', ydata='temp', f1data='conc', f2data='time', indicator='NPV', nsamples=20, figures_path=figures_path)
        analysis_optimization_leaching(system=sys, fununit=fununit, xdata='solventRatio', ydata='time', f1data='conc', f2data='temp', indicator='NPV', nsamples=20, figures_path=figures_path)
        analysis_optimization_leaching(system=sys, fununit=fununit, xdata='solventRatio', ydata='conc', f1data='time', f2data='temp', indicator='NPV', nsamples=20, figures_path=figures_path)
        analysis_optimization_leaching(system=sys, fununit=fununit, xdata='conc', ydata='time', f1data='temp', f2data='solventRatio', indicator='NPV', nsamples=20, figures_path=figures_path)
        analysis_optimization_leaching(system=sys, fununit=fununit, xdata='conc', ydata='temp', f1data='time', f2data='solventRatio', indicator='NPV', nsamples=20, figures_path=figures_path)
        analysis_optimization_leaching(system=sys, fununit=fununit, xdata='temp', ydata='time', f1data='conc', f2data='solventRatio', indicator='NPV', nsamples=20, figures_path=figures_path)

    # Target Analysis
    # --------------------------
    if desire_target == 'yes':
        analysis_target(fununit=fununit, feedPG=feedPG, REEcontent=REEcontent, num_ind_REEs=num_ind_REEs, num_samples=num_samples, desire_target=desire_target, figures_path=figures_path)

    # Scenario Analysis
    # --------------------------
    if desire_scenario == 'yes':
        analysis_scenario(fununit, num_ind_REEs, figures_path)

    return print('-----analysis complete-----')


run_analysis(fununit='PG', feedPG=1000000, REEcontent=0.5/100, num_ind_REEs=9,
             report='no', 
             num_samples=300, uncertainty='no', sensitivity='no', parameter='technological', 
             optimization='no', 
             desire_target='no', 
             desire_scenario='no')

