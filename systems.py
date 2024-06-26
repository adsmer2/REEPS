"""
Rare Earth Element recovery from Phosphogypsum System (REEPS)

The Pennsylvania State University
Chemical Engineering Department
S2D2 Lab (Dr. Rui Shi)
@author: Adam Smerigan
"""
# %%
import collections.abc
#hyper needs the four following aliases to be done manually.
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
collections.MutableSet = collections.abc.MutableSet
collections.MutableMapping = collections.abc.MutableMapping

# import packages
import biosteam as bst
import qsdsan as qs
from qsdsan.utils import clear_lca_registries
import os, pickle, numpy as np, pandas as pd, qsdsan as qs

# import functions from other files
from _component import *
from _lca_data import *
from model import *
from my_TEA import my_TEA

# import unit operations
from leaching import LeachingSulfuric
from filter_gypsum_l import filter_gypsum_l
from filter_gypsum_u import filter_gypsum_u
from precipitation_oxalate import precipitation_oxalate
from filter_oxalate import filter_oxalate
from fired_heater_oxalate import fired_heater_oxalate
from resuspension_tank import resuspension_tank
from black_box import black_box
from precipitation_oxalate_pure import precipitation_oxalate_pure
from filter_oxalate_pure import filter_oxalate_pure
from fired_heater_oxalate_pure import fired_heater_oxalate_pure
from Mixer_wastewater import Mixer_wastewater
from precipitation_NaPhosphate import precipitation_NaPhosphate
from PGmixer import PGmixer
from Fake_mixer import Fake_mixer


bw_path = os.path.dirname(__file__)
data_path = os.path.join(bw_path, 'data')

c_path = os.path.dirname(__file__)
figures_path = os.path.join(c_path, 'figures')

# %%

# =============================================================================
# Create the System
# =============================================================================

# Collect components
cmps = create_components() # generates the properties required for each component used in the simulation
cmps.set_synonym('H2O', 'Water') # set an alternative name for water that can be used

# Get the LCA emission factors and define functional unit
ind_df_processed, cf_dct =_load_lca_data() # Generates impact items that can be linked to flows in/out of the system. Items contain the emission factors calculated from ecoinvent v3.9.1 data for each indicator and flow specified in _lca_data.py
# Function that creates and simulates the system along with the corresponding TEA and LCA results

def create_system(fununit, feedPG, REEcontent, num_ind_REEs):
    """
    fununit = 'PG' # functional unit of interest (PG, REO, or none)

    feedPG = 1000000 #amount of PG being leached (kg/hr)

    REEcontent = 0.5/100 #REE content in PG (wt fraction)

    num_ind_REEs = 9 # Number of REEs recoverable based on a PG REE composition cutoff = 1% for 9 individual REEs valued at 51.5 $/kg REO w/ 97.8% of total REE mass still recoverable. Other values in the comments below
    """
    # --------------------------------------------------------------------------
    # Create inlet/outlet streams and link them to price and environmental impact
    # --------------------------------------------------------------------------
    Nd_item = qs.ImpactItem.get_item('neodymium_item') # For the REO product activity, get the stored environmental impact item generated from _lca_item.py
    PGstack_item = qs.ImpactItem.get_item('PGstack_item') # For the PG stack activity, get the stored environmental impact item generated from _lca_item.py

    # create conditional statements for different scenarios
    # Scenarios for the number of REEs recoverable
    if num_ind_REEs == 14: # number of recoverable REEs based on relative abundance cutoff 0%. Other values in Excel sheet. 2 REEs Lu and Sm were not observed in the sample of PG
        Nd_price = 54.99
        cutoff_mass_REE = 1 # mass fraction of REEs recoverd 
    elif num_ind_REEs == 9: # number of recoverable REEs based on relative abundance cutoff 1%. Other values in Excel sheet 
        Nd_price = 51.5
        cutoff_mass_REE = 0.978 # mass fraction of REEs recoverd 
    elif num_ind_REEs == 7: # number of recoverable REEs based on relative abundance cutoff 2%. Other values in Excel sheet 
        Nd_price = 43.7 # price of REEs recovered
        cutoff_mass_REE = 0.948 # mass fraction of REEs recoverd 
    else:
        raise NameError(f'In systems, {num_ind_REEs} is not a defined value "14", "9", or "7")')
    REEinPG = REEcontent *cutoff_mass_REE # Recoverable REE content in PG (wt fraction)

    # Scenarios for different functional units
    # creates streams with (or without depending on the functional unit) the credits for PG stack treatment and REO production
    if fununit == 'PG': # functional unit of 1 metric ton of PG remediated
        qs.SanStream('Ln2O3', stream_impact_item=Nd_item, price=Nd_price) # this fununit gets a credit for the REO product produced while treating the PG waste
        qs.SanStream('rawPG')
    elif fununit == 'REO': # functional unit of 1 kg of REO produced
        qs.SanStream('Ln2O3', price=Nd_price)
        qs.SanStream('rawPG', stream_impact_item=PGstack_item) # this fununit gets a credit for removing the need for stack treatment of PG while producing REOs
    elif fununit == 'none': # no functional unit (system with no credits)
        qs.SanStream('Ln2O3', price=Nd_price)
        qs.SanStream('rawPG')
    else:
        raise RuntimeError(f'In systems, {fununit} is not "Ln", "PG", or "none"')

    # Gypsum (coproduct)
    gypsum_item = qs.StreamImpactItem.get_item('gypsum_item') # Call the environmental impact item generated from _lca_item.py
    qs.SanStream('gypsum', stream_impact_item=gypsum_item, price=0.00828) # Initialize stream and link emission factors and price
    # Price calculated from US mineral commodities survey 2022

    # Sulfuric Acid
    H2SO4_item = qs.StreamImpactItem.get_item('H2SO4_item') # Call the environmental impact item generated from _lca_item.py
    qs.SanStream('lixiviant_acid', stream_impact_item=H2SO4_item, price=0.089/0.984) # Initialize stream and link emission factors and price. 
    # Sulfuric acid sold typically as 98.4 wt% solution. I assume this stream is 100% sulfuric acid.
    # Ecoinvent models pure substances

    # Oxalic acid
    OA_item = qs.StreamImpactItem.get_item('OA_item') # Call the environmental impact item generated from _lca_item.py
    qs.SanStream('OA_feed', stream_impact_item=OA_item, price=0.855) # Initialize stream and link emission factors and price
    OA_item_2 = qs.StreamImpactItem.get_item('OA_item_2') # Call the environmental impact item generated from _lca_item.py
    qs.SanStream('OA_feed_2', stream_impact_item=OA_item_2, price=0.855) # Initialize stream and link emission factors and price

    # Process Water
    pwater_item = qs.StreamImpactItem.get_item('pwater_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.SanStream('lixiviant_water', stream_impact_item=pwater_item, price=0.000367) # Initialize stream and link emission factors and price
    pwater_item_2 = qs.StreamImpactItem.get_item('pwater_item_2') # Call the environmental impact item generated from _lca_item.py
    qs.SanStream('RS_water', stream_impact_item=pwater_item_2, price=0.000367) # Initialize stream and link emission factors and price
    # price from pg 64. Product and Process Design Principles​. Synthesis, Analysis and Design, Third Edition - Seider, Seader, Lewin, Widagdo

    electricity_item = qs.StreamImpactItem.get_item('electrictiy_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.PowerUtility.price = 0.0976 # $/(kWh) electricity price for industry in Florida in July 2022 https://www.eia.gov/electricity/data/browser/#/topic/7?agg=0,1&geo=0000001&endsec=2&linechart=ELEC.PRICE.FL-IND.M&columnchart=ELEC.PRICE.FL-IND.M&map=ELEC.PRICE.FL-IND.M&freq=M&ctype=columnchart&ltype=pin&rtype=s&maptype=0&rse=0&pin=
    # Used in the qs.LCA() input below to quantify impacts as the system changes
    
    heatNG_item = qs.StreamImpactItem.get_item('heatNG_item') # Call the environmental impact item generated from _lca_item.pyp
    # Used in the qs.LCA() input below to quantify impacts as the system changes
    # price_NG = 7.534*10**(-6) # $/kJ. Price of NG. NG price for industry July 2022 from EIA
    # The cost of using natural gas is accounted for in the individual unit that uses it by including self.add_OPEX = {'Natural Gas': price_NG*self.heat_duty}.
    

    CO2_item = qs.StreamImpactItem.get_item('CO2_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.SanStream('CO2_emission', stream_impact_item=CO2_item) # Initialize stream and link emission factors 
    CO2_item_2 = qs.StreamImpactItem.get_item('CO2_item_2') # Call the environmental impact item generated from _lca_item.pyp
    qs.SanStream('CO2_emission_2', stream_impact_item=CO2_item_2) # Initialize stream and link emission factors

    CO_item = qs.StreamImpactItem.get_item('CO_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.SanStream('CO_emission', stream_impact_item=CO_item) # Initialize stream and link emission factors
    CO_item_2 = qs.StreamImpactItem.get_item('CO_item_2') # Call the environmental impact item generated from _lca_item.pyp
    qs.SanStream('CO_emission_2', stream_impact_item=CO_item_2) # Initialize stream and link emission factors

    wastewater_item = qs.StreamImpactItem.get_item('wastewater_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.WasteStream('wastewater_treatment', stream_impact_item=wastewater_item) # Initialize stream and link emission factors

    HNO3_item = qs.StreamImpactItem.get_item('HNO3_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.WasteStream('HNO3_feed', stream_impact_item=HNO3_item, price=0.393/0.98) # Initialize stream and link emission factors and price
    # Nitric acid typically sold as 98% solution. I assume this stream is 100% nitric acid.
    # Ecoinvent models pure substances

    NaOH_item = qs.StreamImpactItem.get_item('NaOH_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.WasteStream('NaOH_feed', stream_impact_item=NaOH_item, price=0.380) # Sold as caustic soda flakes
    # Ecoinvent models pure substances

    Na3PO4_item = qs.StreamImpactItem.get_item('Na3PO4_item') # Call the environmental impact item generated from _lca_item.pyp
    qs.WasteStream('Na3PO4_feed', stream_impact_item=Na3PO4_item, price=0.505) # Initialize stream and link emission factors and price

    # --------------------------------------------------------------------------
    # Create the System
    # --------------------------------------------------------------------------

    #use qs.Flowsheet to connect streams linked to LCA impact items to actual streams in the system
    flowsheet = qs.Flowsheet.flowsheet.default # shortcut to call the system flowsheet
    fs_stream = flowsheet.stream # shortcut to call specific streams
    fs_unit = flowsheet.unit # shortcut to call specific unit operations
    
    # Set up each unit in the system
    # must initialize with a unit ID, inlet streams, outlet streams, and whatever other attributes are coded into the unit python script
    M0 = PGmixer(ID='M0', ins=fs_stream.rawPG, outs='PG', feedPG=feedPG, REEinPG=REEinPG, UinPG=0.003159/100) # dummy unit (mixer) to make data interpretation easier later on
    U1 = LeachingSulfuric(ins=(M0-0, fs_stream.lixiviant_water, fs_stream.lixiviant_acid, 'rec_lixiviant'), outs=('leachate', 'underflow'), feedPG=feedPG, REEinPG=REEinPG) # leaching
    F1 = filter_gypsum_l(ID='F1', ins=(U1.outs[0]), outs=('solids', 'dilute_Ln')) # filter (overflow)
    F2 = filter_gypsum_u(ID='F2', ins=(U1.outs[1]), outs=('u_solids', 3-U1), Ln_desired=False) # filter (underflow)
    P1 = precipitation_oxalate(ID='P1', ins=(F1-1, fs_stream.OA_feed), outs=('conc_Ln_oxalate', 'wastewater','Hion_1')) # oxalate precipitation
    F3 = filter_oxalate(ID='F3', ins=P1-0, outs=('Ln_oxalate', 'wastewater_2')) # filter (mixed REE-oxalate)
    P3 = precipitation_NaPhosphate(ID='P3', ins=(P1-1, fs_stream.Na3PO4_feed, fs_stream.NaOH_feed, F3-1), outs=('U_conc', 'wastewater_5')) # radionuclide precipitation 
    H1 = fired_heater_oxalate(ID='H1', ins=F3-0, outs=('mixed_Ln2O3_solid', fs_stream.CO2_emission, fs_stream.CO_emission)) # fired heater (oxalate removal)
    RS = resuspension_tank(ID= 'RS', ins=(H1-0,fs_stream.HNO3_feed, fs_stream.RS_water), outs=('conc_Ln')) # Mixed REO redissolving tank
    S1 = black_box(ID= 'S1', ins=(RS-0), outs=('pure_Ln', 'Ln_lost')) # Selective REE separation
    P2 = precipitation_oxalate_pure(ID='P2', ins=(S1-0, fs_stream.OA_feed_2), outs=('pure_Ln_oxalate', 'wastewater_3','Hion_2'), num_parallel_units=num_ind_REEs) # oxalate precipitation
    F4 = filter_oxalate_pure(ID='F4', ins=(P2-0), outs=('pure_Ln_oxalate_s', 'wastewater_4'), num_parallel_units=num_ind_REEs) # filter (individual REE-oxalates)
    H2 = fired_heater_oxalate_pure(ID='H2', ins=F4-0, outs=('pure_Ln2O3', fs_stream.CO2_emission_2, fs_stream.CO_emission_2)) # # fired heater (oxalate removal and calcination)
    M1 = Fake_mixer(ID='M1', ins=(F1-0, F2-0), outs=fs_stream.gypsum) # dummy unit (mixer) to make data interpretation easier later on
    M4 = Fake_mixer(ID='M4', ins=H2-0, outs=fs_stream.Ln2O3) # dummy unit (mixer) to make data interpretation easier later on
    WT = Mixer_wastewater(ID='WT', ins=(P3-1, P2-1, F4-1), outs=(fs_stream.wastewater_treatment), init_with='WasteStream') # wastewater treatment plant
    
    # qs.sanunits.Mixer
    # bst.units.Mixer
    # simulate each unit in the system
    M0.simulate()
    U1.simulate()
    F1.simulate()
    F2.simulate()
    P1.simulate()
    P3.simulate()
    F3.simulate()
    H1.simulate()
    RS.simulate()
    S1.simulate()
    P2.simulate()
    F4.simulate()
    H2.simulate()
    M1.simulate()
    M4.simulate()
    WT.simulate()

    # Initialize the system
    sys = qs.System('sys', path=(M0,U1,F1,F2,M1,P1,P3,F3,H1,RS,S1,P2,F4,H2,M4,WT), recycle=F2-1) # initialize system and include recycle stream for leaching lixiviant
       
    # -------------------------------------------------------------------------
    # Calculate TEA and LCA
    # -------------------------------------------------------------------------

    # Labor Costs
    # Calculated according to Seider - 4th edition product and process design principles
    operators_per_section = 2 # operators per section from Seider recommendation
    num_process_sections = 2 # number of proces sections from Seider recommendation
    num_operators_per_shift = operators_per_section * num_process_sections *2 # multiplied by 2 for large continuous flow process (e.g., 1000 ton/day product). from Seider pg 505
    num_shifts = 5 # number of shifts 
    pay_rate = 50 # $/hr
    DWandB = num_operators_per_shift*num_shifts*2080*pay_rate # direct wages and benefits. DWandB [$/year] = (operators/shift)*(5 shifts)*(40 hr/week)*(operating days/year-operator)*($/hr)
    Dsalaries_benefits = 0.15*DWandB # direct salaries and benefits from Seider
    O_supplies = 0.06*DWandB # Operating supplies and services from Seider
    technical_assistance = 5*75000 # $/year. Technical assistance to manufacturing. assume 5 workers at $75000/year 
    control_lab = 5*80000 # $/year. Control laboratory. assume 5 workers at $80000/year 
    labor = DWandB + Dsalaries_benefits + O_supplies + technical_assistance + control_lab

    # Using CEPCI for June 2022 of 832.6 (Seider et al from 2006 has CEPCI of 500) for a CEPCI adjustment ratio of 1.6652
    # All equipment purchase prices from the Seider 2006 book are adjusted accordingly (multiplied by 1.6652)
    tea = my_TEA(system=sys, # Should contain all feed and product streams
                IRR=0.15, # Interest rate (fraction) Seider 4th ed page 526 
                duration=(2023, 2053), # Start and end year of venture (e.g. (2018, 2038))
                depreciation='MACRS7', # 'MACRS' + number of years (e.g. 'MACRS7'). Using 7 years based on IRS Asset Class 49.5, “Waste Reduction and Resource Recovery Plants.”. NREL Davis - 2016
                income_tax=0.28, # Combined federal and state income tax rate (fraction). 21% Fed and assume an avg state tax rate ~7%
                operating_days=328, # Number of operating days per year. Assumed 90%
                lang_factor=4.28, # Lang factor for getting fixed capital investment from total purchase cost. If no lang factor, estimate capital investment using bare module factors. From Seider - 4th edition pg 447 for solids-fluids processing plant
                construction_schedule=(0.08, 0.60, 0.32), # Startup investment fractions per year (e.g. (0.5, 0.5) for 50% capital investment in the first year and 50% investment in the second). From NREL TEA reports
                WC_over_FCI=0.05, # Working capital as a fraction of fixed capital investment. 5% used in NREL TEA reports
                maintenance=0.045, # Yearly fee as a fraction of total depreciable capital. From Seider - 4th edition pg 500
                labor = labor, # calculated above using similar method to Seider 4th ed
                membrane_cost=fs_unit.S1.membrane_cost # $. The purchase cost of the membrane adsorption unit. This attribute bypasses the lang factor for this expense.
   )
    
    # lca calculates electricity and heat impacts separate from stream impacts. 
    # Electricity has fununit of kW*h and normalized to lifetime of a year. Automatically updates with additions to model
    # Heat from NG has fununit of MJ and is normalized to the lifetime. Must be manually updated with additions to model
    lca = qs.LCA(system=sys, lifetime=1, lifetime_unit='year', 
                 uptime_ratio = tea.operating_days/365, 
                 electricity_item=lambda: sys.power_utility.rate*(tea.operating_days*24), 
                 heatNG_item=lambda:(U1.heat_duty + H1.heat_duty + H2.heat_duty)/1000*(tea.operating_days*24))
    
    return sys, lca, tea
