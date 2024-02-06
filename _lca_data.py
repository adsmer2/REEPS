"""
Phosphogypsum Rare Earth Element Recovery System (PREERS)

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""

# %%

import os, sys, pickle, pandas as pd, qsdsan as qs, brightway2 as bw
import re as re

c_path = os.path.dirname(__file__)
data_path = os.path.join(c_path, 'data')

__all__ = ('create_indicators', 'create_items', 'get_cf_data', 'save_cf_data', '_load_lca_data')

# %%
# =======================================================================
# Create/open a brightway2 project
# =======================================================================
bw.projects.set_current('REEPSv1') # create/open brightway2 project
bw.bw2setup() # setup the project and make sure it has the biosphere3 database
bsdb = bw.Database('biosphere3') # biosphere database for ecoinvent 3.9.1


# =============================================================================
# Database download, only run this once to download the database
# =============================================================================

# To use this function, you must have purchased access to ecoinvent v3.9.1
# If this doesn't work, try running eidl in jupyter nb 
def download_data():
    import eidl
    eidl.get_ecoinvent() # ei cutoff 3.9.1

# download_data() # uncomment this to import the ecoinvent database. Comment it again after importing
eidb = bw.Database('cutoff391') # LCI database from ecoinvent (v3.9.1 cutoff method). Contains technological flows and activities

# =============================================================================
# Using elementary flows in the biosphere3 database
# =============================================================================
def create_biosphere_activity():
    # In order to use Biosphere3 database flows, they must be added to the Ecoinvent database via the creation of a new activity with one exchange (the bsdb flow of interest)
    # Loading activities to bw2qsd in "select_items()" and then calculating emission factors in "organize_cfs()" is good after doing the above step.

    for activity in [act for act in eidb if act['name']=='CO2_emission' or act['name']=='CO_emission' or act['name']=='water_emission']:
        activity.delete()

    # create new activity in the eidb with the name of your biosphere flow
    CO2_emission = eidb.new_activity(code = 'One CO2 stream=1', name = "CO2_emission", unit = "kg")
    CO2_emission.save()

    # select the biosphere flow of interest and add it as a new exchange to the activity created
    CO2_fossil = [i for i in bsdb if 'Carbon dioxide, fossil' in i['name'] and 'low population density, long-term' in i['categories']][0]
    CO2_emission.new_exchange(input=CO2_fossil.key,amount=1,type='biosphere').save()

    # confirm your flow is now in the eidb database correctly
    confirm_CO2 = [i for i in eidb if 'CO2_emission' in i['name']][0]
    # -----------
    # create new activity in the eidb with the name of your biosphere flow
    CO_emission = eidb.new_activity(code = 'One CO stream=1', name = "CO_emission", unit = "kg")
    CO_emission.save()

    # select the biosphere flow of interest and add it as a new exchange to the activity created
    CO_fossil = [i for i in bsdb if 'Carbon monoxide, fossil' in i['name'] and 'low population density, long-term' in i['categories']][0]
    CO_emission.new_exchange(input=CO_fossil.key,amount=1,type='biosphere').save()

    # confirm your flow is now in the eidb database correctly
    confirm_CO = [i for i in eidb if 'CO_emission' in i['name']][0]

    # # -----------
    # # create new activity in the eidb with the name of your biosphere flow
    # water_emission = eidb.new_activity(code = 'One water_surface stream=1', name = "water_emission", unit = "m3")
    # water_emission.save()

    # # select the biosphere flow of interest and add it as a new exchange to the activity created
    # water_surface = [i for i in bsdb if 'Water' in i['name'] and 'water' in i['categories'] and 'surface water' in i['categories'] and 'emission' in i['type']][0]
    # water_emission.new_exchange(input=water_surface.key,amount=1,type='biosphere').save()

    # # confirm your flow is now in the eidb database correctly
    # confirm_water = [i for i in eidb if 'water_emission' in i['name']][0]
    return confirm_CO2, confirm_CO #, confirm_water

# ==========================================================================
# Adding activity into ecoinvent database for PG stacking
# ==========================================================================
def create_PGstack_activity():
    # Loading custom activities to bw2qsd in "select_items()" and then calculating emission factors in "organize_cfs()" is good after doing the above step.

    # create new activity in the eidb with the name of your custom activity
    PGstack_Tsioka = eidb.new_activity(code = 'PG Stack Tsioka and Voudrias', name = "PG Stack Tsioka", unit = "t") # t is tonne in bw2. ('t', 'kilogram', 1e3) from https://github.com/brightway-lca/brightway2-io/blob/main/bw2io/units.py
    for activity in [act for act in eidb if act['name']=='PG Stack Tsioka']:
        activity.delete()
    PGstack_Tsioka.save()
    # Add exchanges to the activity according to the paper, Journal of Cleaner Production 266 (2020) 121386, Tsioka and Voudrias - 2020 - Comparison of alternative management methods for phosphogypsum waste using life cycle analysis
    land_occupation = [exc for exc in bsdb if 'Occupation, pasture, man made' in exc['name'] and not 'intensive' in exc['name'] and not 'extensive' in exc['name']][0]
    PGstack_Tsioka.new_exchange(input=land_occupation.key,amount=153.8,type='biosphere').save()
    fluoride = [exc for exc in bsdb if 'Fluoride' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=fluoride.key,amount=118.12/1000,type='biosphere').save() 
    chromium = [exc for exc in bsdb if 'Chromium VI' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0] # chromium IV gives an ecotoxicity closer to that of the Tsioka and Voudrias paper (CrIII=0Pt, CrIV=0.0046Pt to literature=0.006Pt)
    PGstack_Tsioka.new_exchange(input=chromium.key,amount=0.000862,type='biosphere').save()
    zinc = [exc for exc in bsdb if 'Zinc II' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=zinc.key,amount=0.00384,type='biosphere').save()
    cadmium = [exc for exc in bsdb if 'Cadmium II' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=cadmium.key,amount=0.00189,type='biosphere').save()
    copper = [exc for exc in bsdb if 'Copper' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=copper.key,amount=0.000794,type='biosphere').save()
    phosphate = [exc for exc in bsdb if 'Phosphate' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=phosphate.key,amount=0.726,type='biosphere').save()
    radium = [exc for exc in bsdb if 'Radium-226' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=radium.key,amount=9.27,type='biosphere').save()
    sulfate = [exc for exc in bsdb if 'Sulfate' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=sulfate.key,amount=64.59,type='biosphere').save()
    calcium = [exc for exc in bsdb if 'Calcium II' in exc['name'] and 'water' in exc['categories'] and 'ground-' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=calcium.key,amount=12.08,type='biosphere').save()
    hydrogen_fluoride = [exc for exc in bsdb if 'Hydrogen fluoride' in exc['name'] and 'air' in exc['categories'] and 'low population density, long-term' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=hydrogen_fluoride.key,amount=38.4/1000,type='biosphere').save() 
    silicon_tetrafluoride = [exc for exc in bsdb if 'Silicon tetrafluoride' in exc['name'] and 'air' in exc['categories'] and 'low population density, long-term' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=silicon_tetrafluoride.key,amount=49.9/1000,type='biosphere').save() # ,unit='gram'
    particulates_10um = [exc for exc in bsdb if 'Particulate Matter, > 10' in exc['name'] and 'air' in exc['categories'] and 'low population density, long-term' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=particulates_10um.key,amount=0.696/1000,type='biosphere').save() 
    uranium_air = [exc for exc in bsdb if 'Uranium-238' in exc['name'] and 'air' in exc['categories'] and 'low population density, long-term' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=uranium_air.key,amount=0.024/1000,type='biosphere').save() 
    radium_air = [exc for exc in bsdb if 'Radium-226' in exc['name'] and 'air' in exc['categories'] and 'low population density, long-term' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=radium_air.key,amount=0.358/1000,type='biosphere').save() 
    radon_air = [exc for exc in bsdb if 'Radon-222' in exc['name'] and 'air' in exc['categories'] and 'low population density, long-term' in exc['categories']][0]
    PGstack_Tsioka.new_exchange(input=radon_air.key,amount=3.92060*10**5,unit='kBq',type='biosphere').save()

    # confirm that the activity saved to the eidb database
    confirm1 = [i for i in eidb if 'PG Stack Tsioka' in i['name']]
    # confirm that exchanges saved to the activity correctly
    confirm2 = [i for i in PGstack_Tsioka.exchanges()]
    return confirm1, confirm2
    
# %%

# =============================================================================
# Impact indicators
# =============================================================================

def create_indicators(replace=True):
    from bw2qsd import CFgetter
    from bw2qsd.utils import format_name
    cutoff391 = CFgetter('cutoff391')
    # ecoinvent version 3.9.1, cutoff
    cutoff391.load_database('cutoff391')

    # Include ReCiPe 2016 v1.03 (Hierarchist) Midpoint LCIA method w/ LT impact
    cutoff391.load_indicators(add=True, method='ReCiPe 2016 v1.03, midpoint (H)', method_exclude=('no LT'))

    # Make the names of the indicators nicer
    ind_df_raw = cutoff391.export_indicators(show=False, path='') # writes indicators to file
    ind_df_processed = ind_df_raw.copy()

    replace = True
    ind_df_processed = ind_df_raw.copy()

    for num, ind in ind_df_raw.iterrows():
        old_name = ind['indicator']
        if 'ReCiPe' in ind['method']: # need to remove special characters for .load_from_file to work
            name1 = old_name.split('(')[1]
            new_name = name1.split(')')[0]
            
        else:
            RuntimeError('need to write code for indicator names of other LCIA methods')
        ind_df_processed.iloc[num]['indicator'] = new_name

    ind_df_processed.sort_values(by=['method', 'category', 'indicator'], inplace=True)

    # overwrites impact indicators if requested
    if replace:
        for ind_ID in ind_df_processed.indicator:
            ind = qs.ImpactIndicator.get_indicator(ind_ID)
            if ind:
                stdout = sys.stdout
                sys.stdout = open(os.devnull, 'w')
                ind.deregister()
                sys.stdout = stdout

    qs.ImpactIndicator.load_from_file(ind_df_processed) # get the stored impact assessment method impact indicators from file
    indicators = qs.ImpactIndicator.get_all_indicators() # link these indicators to the qsdsan package for use in the system model

    return cutoff391, ind_df_processed, indicators


# %%

# =============================================================================
# Impact items
# =============================================================================

def select_items(database):
    '''
    Identify the activities of interest from ecoinvent using brightway2
    '''
    all_acts = {}
    eidb = bw.Database('cutoff391')

    def new_act(database, all_acts, name):
        act = database.copy(name)
        all_acts[name] = act
        return act

    # Catch all printouts (i.e., don't show them in the console)
    stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    # Sulfuric Acid
    H2SO4_item = new_act(database, all_acts, 'H2SO4_item')
    # H2SO4_item.load_activities('market for sulfuric acid', add=True, filter={'location': 'RoW'}, limit=None)
    # to_remove = []
    # for act in H2SO4_item.activities.keys():
    #     if ('bromine production' in act) or ('sodium chloride production, powder' in act):
    #         to_remove.append(act)
    # H2SO4_item.remove('activity', to_remove)
    dct = [i for i in eidb if 'market for sulfuric acid' in i['name'] and 'RoW' in i['location']][0]
    dct2 = {dct['name']: dct}
    H2SO4_item._activities.update(dct2)

    # Oxalic Acid
    OA_item = new_act(database, all_acts, 'OA_item')
    dct = [i for i in eidb if 'market for oxalic acid' in i['name'] and 'GLO' in i['location']][0]
    dct2 = {dct['name']: dct}
    OA_item._activities.update(dct2)

    # Process Water
    pwater_item = new_act(database, all_acts, 'pwater_item')
    dct = [i for i in eidb if 'market for water, decarbonised' in i['name'] and 'US' in i['location']][0]
    dct2 = {dct['name']: dct}
    pwater_item._activities.update(dct2)

    # Generic Wastewater Treatment (represents a negative input, treatment of a waste)
    wastewater_item = new_act(database, all_acts, 'wastewater_item')
    dct = [i for i in eidb if 'market for wastewater, average' in i['name'] and 'RoW' in i['location']][0]
    dct2 = {dct['name']: dct}
    wastewater_item._activities.update(dct2)

    # Gypsum
    gypsum_item = new_act(database, all_acts, 'gypsum_item')
    dct = [i for i in eidb if 'market for gypsum, mineral' in i['name'] and 'RoW' in i['location']][0]
    dct2 = {dct['name']: dct}
    gypsum_item._activities.update(dct2)

    # Neodymium
    neodymium_item = new_act(database, all_acts, 'neodymium_item')
    dct = [i for i in eidb if 'market for neodymium' in i['name'] and 'GLO' in i['location']][0]
    dct2 = {dct['name']: dct}
    neodymium_item._activities.update(dct2)

    # rare earth oxides production, from rare earth oxide concentrate, 70% REO
    REE_prod_item = new_act(database, all_acts, 'REE_prod_item')
    dct = [i for i in eidb if 'rare earth oxides production, from rare earth oxide concentrate, 70% REO' in i['name']][0]
    dct2 = {dct['name']: dct}
    REE_prod_item._activities.update(dct2)

    # Electricity
    electricity_item = new_act(database, all_acts, 'electricity_item')
    dct = [i for i in eidb if 'market for electricity, medium voltage' in i['name'] and 'US-SERC' in i['location']][0]
    dct2 = {dct['name']: dct}
    electricity_item._activities.update(dct2)

    # Steel
    steel_item = new_act(database, all_acts, 'steel_item')
    dct = [i for i in eidb if 'market for steel, chromium steel 18/8'==i['name']][0]
    dct2 = {dct['name']: dct}
    steel_item._activities.update(dct2)

    # Heat Production, Natural Gas
    heatNG_item = new_act(database, all_acts, 'heatNG_item')
    dct = [i for i in eidb if 'heat production, natural gas, at industrial furnace >100kW' in i['name'] and 'CA-QC' in i['location']][0]
    dct2 = {dct['name']: dct}
    heatNG_item._activities.update(dct2)

    # Carbon Dioxide
    CO2_item = new_act(database, all_acts, 'CO2_item')
    dct = [i for i in eidb if 'CO2_emission' in i['name']][0]
    dct2 = {dct['name']: dct}
    CO2_item._activities.update(dct2)

    # Carbon Monoxide
    CO_item = new_act(database, all_acts, 'CO_item')
    dct = [i for i in eidb if 'CO_emission' in i['name']][0]
    dct2 = {dct['name']: dct}
    CO_item._activities.update(dct2)

    # # Wastwater treatment: medium density fibreboard
    # wastewater_pb_item = new_act(database, all_acts, 'wastewater_pb_item')
    # dct = [i for i in eidb if 'treatment of wastewater from medium density fibreboard production' in i['name'] and 'RER' in i['location']][0]
    # dct2 = {dct['name']: dct}
    # wastewater_pb_item._activities.update(dct2)

    # Nitric Acid
    HNO3_item = new_act(database, all_acts, 'HNO3_item')
    dct = [i for i in eidb if 'market for nitric acid' in i['name'] and 'RoW' in i['location']][0]
    dct2 = {dct['name']: dct}
    HNO3_item._activities.update(dct2)

    # Sodium Hydroxide
    NaOH_item = new_act(database, all_acts, 'NaOH_item')
    dct = [i for i in eidb if 'market for sodium hydroxide' in i['name'] and 'GLO' in i['location']][0] # sodium hydroxide, without water, in 50% solution state
    dct2 = {dct['name']: dct}
    NaOH_item._activities.update(dct2)

    # Sodium Phosphate
    Na3PO4_item = new_act(database, all_acts, 'Na3PO4_item')
    dct = [i for i in eidb if 'market for trisodium phosphate' in i['name'] and 'GLO' in i['location']][0] # anhydrous trisodium phosphate
    dct2 = {dct['name']: dct}
    Na3PO4_item._activities.update(dct2)

    # PG stack
    PGstack_item = new_act(database, all_acts, 'PGstack_item')
    dct = [i for i in eidb if 'PG Stack Tsioka' in i['name']][0] # Custom PGstack activity
    dct2 = {dct['name']: dct}
    PGstack_item._activities.update(dct2)

    # # Water, Emission to surface
    # water_emission_item = new_act(database, all_acts, 'water_emission_item')
    # dct = [i for i in eidb if 'water_emission' in i['name']][0] # Custom eidb water_emission activity from flow from biosphere3 database
    # dct2 = {dct['name']: dct}
    # water_emission_item._activities.update(dct2)

    # Restore printouts
    sys.stdout = stdout

    return all_acts 


def get_stats(df, keep_raw_data=False):
    '''
    makes the raw impact data more organized
    '''
    df2 = pd.DataFrame(columns=df.columns, dtype='float64') if not keep_raw_data else df.copy()
    # df2 = df2.append(df[1:].min(), ignore_index=True)
    # df2 = df2.append(df[1:].mean(), ignore_index=True)
    # df2 = df2.append(df[1:].max(), ignore_index=True)
    df2 = pd.concat([df[1:].min(), df[1:].apply(pd.to_numeric, errors='coerce').mean(), df[1:].max()], axis=1, ignore_index=True).transpose()
    df2.loc[-3:, ('-', '-', 'activity name')] = ('min', 'mean', 'max')
    functional_unit = df.loc[1, ('-', '-', 'functional unit')]
    df2.loc[1:, ('-', '-', 'functional unit')] = functional_unit
    return df2

def organize_cfs(all_acts):
    '''
    Apply get_stats to the data and make any changes to the emission factors necesary before writing to file
    '''
    cf_dct = {}
    cf_dct['H2SO4_item'] = get_stats(all_acts['H2SO4_item'].CFs)

    cf_dct['OA_item'] = get_stats(all_acts['OA_item'].CFs)
    cf_dct['OA_item_2'] = get_stats(all_acts['OA_item'].CFs)

    cf_dct['steel_item'] = get_stats(all_acts['steel_item'].CFs)

    cf_dct['pwater_item'] = get_stats(all_acts['pwater_item'].CFs)
    cf_dct['pwater_item_2'] = get_stats(all_acts['pwater_item'].CFs)

    cf_dct['gypsum_item'] = get_stats(all_acts['gypsum_item'].CFs)
    cols = all_acts['gypsum_item'].CFs.columns[2:]
    cf_dct['gypsum_item'][cols] *= -1 # credit

    cf_dct['wastewater_item'] = get_stats(all_acts['wastewater_item'].CFs)
    cols = all_acts['wastewater_item'].CFs.columns[2:]
    cf_dct['wastewater_item'][cols] *= -1/997 # turning negative input into output. Divide CF (impact/m3 wastewater) by 997 to change functional unit from 1 m3 to 1 kg, since 997 kg are in 1 m3. Note: only changes in the values of the CF reflected in Excel (not the FU definition in text)

    cf_dct['neodymium_item'] = get_stats(all_acts['neodymium_item'].CFs)
    cols = all_acts['neodymium_item'].CFs.columns[2:]
    cf_dct['neodymium_item'][cols] *= -1 # credit

    cf_dct['REE_prod_item'] = get_stats(all_acts['REE_prod_item'].CFs)
    cols = all_acts['REE_prod_item'].CFs.columns[2:]
    cf_dct['REE_prod_item'][cols] *= -1 # credit

    cf_dct['electricity_item'] = get_stats(all_acts['electricity_item'].CFs)

    cf_dct['heatNG_item'] = get_stats(all_acts['heatNG_item'].CFs)

    cf_dct['CO2_item'] = get_stats(all_acts['CO2_item'].CFs)

    all_acts['CO_item'].CFs.loc[1:,('ReCiPe 2016 v1.03, midpoint (H)','climate change','global warming potential (GWP1000)')] = 1.9 # change the GWP100 of carbon monoxide from 0 to 1.9 to account for the indirect GW affect of carbon monoxide from https://archive.ipcc.ch/publications_and_data/ar4/wg1/en/ch2s2-10-3-2.html
    cf_dct['CO_item'] = get_stats(all_acts['CO_item'].CFs)

    # cf_dct['wastewater_pb_item'] = get_stats(all_acts['wastewater_pb_item'].CFs)
    # cols = all_acts['wastewater_pb_item'].CFs.columns[2:]
    # cf_dct['wastewater_pb_item'][cols] *= -1 # turning negative input into output

    cf_dct['HNO3_item'] = get_stats(all_acts['HNO3_item'].CFs)

    cf_dct['NaOH_item'] = get_stats(all_acts['NaOH_item'].CFs) 

    cf_dct['Na3PO4_item'] = get_stats(all_acts['Na3PO4_item'].CFs)

    cf_dct['PGstack_item'] = get_stats(all_acts['PGstack_item'].CFs)
    cols = all_acts['PGstack_item'].CFs.columns[2:]
    cf_dct['PGstack_item'][cols] *= -1/1000 # Turn emission factors negative to represent a waste treatment activity. Change functional unit from metric tonnes to kg. Note: fununit change not reflected in Excel

    # all_acts['water_emission_item'].CFs.loc[1:,('ReCiPe Midpoint (H) V1.13','water depletion','WDP')] = -1 # change the WDP of water emission from 0 to -1 (kg water-eq/kg water) to account for water returning to the environnment

    return cf_dct

def create_items(ind_df_processed, cf_dct, replace=True):
    '''
    creates impact items in qsdsan from the organized emission factors in cf_dct from organize_cfs()
    '''
    items = []
    for item_ID, df in cf_dct.items():
        item = qs.ImpactItem.get_item(item_ID)
        if not (replace and item):
            if not ('item' in item_ID and item_ID!='electricity_item'):
                item = qs.ImpactItem(ID=item_ID,
                                     functional_unit=df.loc[1, ('-', '-', 'functional unit')])
            else:
                item = qs.StreamImpactItem(ID=item_ID)

        for num in ind_df_processed.index:
            ind_ID = ind_df_processed.iloc[num]['indicator']
            ind_str = ind_df_processed.iloc[num]['full_name']
            ind_col = tuple(i for i in ind_str.split("'") if len(i)>2)
            item.add_indicator(ind_ID, CF_value=df[df.values=='mean'][ind_col].item())

        items.append(item)
    return items


    # %%

# =============================================================================
# Run and save data
# =============================================================================

def get_cf_data():
    '''
    Runs all the above functions to prep the data for writing to file in save_cf_data()
    '''
    from bw2qsd import remove_setups_pickle
    cutoff391, ind_df_processed, indicators = create_indicators()
    all_acts = select_items(cutoff391)
    cf_dct = organize_cfs(all_acts)

    # Only run this at the very end to remove the outdated setup.pickle file
    remove_setups_pickle()

    return ind_df_processed, all_acts, cf_dct


def save_cf_data():
    '''
    For each activity in select_items(), writes the emission factors for each indicator
    '''
    ind_df_processed, all_acts, cf_dct = get_cf_data()

    ind_file = 'indicators.tsv'
    raw_path = os.path.join(data_path, 'CFs')

    # Impact indicators
    ind_df_processed.to_csv(os.path.join(data_path, ind_file), sep='\t')

    # Raw data
    if not os.path.isdir(raw_path):
        os.mkdir(raw_path)
    for k, v in all_acts.items():
        v.CFs.to_csv(os.path.join(raw_path, f'CFs_{k}.tsv'), sep='\t')

    # Organized data
    f = open(os.path.join(data_path, 'cf_dct.pckl'), 'wb')
    pickle.dump(cf_dct, f)
    f.close()

# =============================================================================
# Collect Saved LCA Data 
# =============================================================================

def _load_lca_data():
    '''
    Load impact indicator and impact item data from file and create impact items in qsdsan for use in systems.py
    '''
    indicator_path = os.path.join(data_path, 'indicators.tsv')
    indel_col = 0
    ind_df_processed = pd.read_csv(indicator_path, sep='\t', index_col=indel_col)
    qs.ImpactIndicator.load_from_file(indicator_path)

    item_path = os.path.join(data_path, 'cf_dct.pckl')
    f = open(item_path, 'rb')
    cf_dct = pickle.load(f)
    f.close()
    create_items(ind_df_processed, cf_dct)


# confirm_CO2, confirm_CO = create_biosphere_activity()
# print(confirm_CO2, confirm_CO)
# confirm1, confirm2 = create_PGstack_activity() # run ONLY the first time you run _LCA_data.py
# print(confirm1, confirm2)

# save_cf_data() # run the first time and every time you make a change to the LCA items

# Errors (and solutions)
# Pickle deserialization error in file 'C:\Users\ajs8911\AppData\Local\pylca\Brightway3\default.c21f969b5f03d33d43e04f8f136e7682\setups.pickle'
#       Go to the the file destination and delete the setups.pickle file
# OMP: Error #15: Initializing libiomp5md.dll, but found libiomp5md.dll already initialized.
#       manually deleted the libiomp5.dll file
#       file is in your anaconda3/envs/XXX/Library/bin folder where XXX is your env name or just anaconda3/Library/bin if not in an env
#       fix from https://stackoverflow.com/questions/20554074/sklearn-omp-error-15-initializing-libiomp5md-dll-but-found-mk2iomp5md-dll-a
