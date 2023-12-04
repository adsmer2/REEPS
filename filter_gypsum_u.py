"""
Phosphogypsum Rare Earth Element Recovery System (PREERS)

Last Modified: 11/04/22

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""
# %%
import qsdsan as qs
import numpy as np

# Create the class
class filter_gypsum_u(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    filt_rate = 1220.5, # filtration rate for coarse particles in kg/hr/m2 from Seider pg 585
    vacuum_energy = 12, # kW/m2. Pressure required to create vacuum (from 2-15 kW/m2) from Hendriksson, Brandt - 2000 - focus on separation in the mining industry
    air_flow = 1.5, # normal m3/m2 filter area/min. air flow through filter cake for pump to withdraw (0.3-2.0 m3/m2 filter area/min) from Hendriksson, Brandt - 2000 - focus on separation in the mining industry
    Ln_desired = True # True if Ln is desired in the liquid stream
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes
        self.filt_rate = filt_rate
        self.vacuum_energy = vacuum_energy
        self.Ln_desired = Ln_desired
        self.air_flow = air_flow

        self.F_BM = {'Rotary Drum Filter': 2.5} # need to confirm this
        
    _N_ins = 1
    _N_outs = 2

    def _run(self):
            
        # Calculate filtration size factor (filter area) 
        self.inlet_solids = self.ins[0].imass['Gypsum'] # kg/hr. Mass of solid gypsum entering the unit
        self.filter_area = self.inlet_solids/self.filt_rate # m2. filter area required

        # Define Outlet Flows
        self.outs[0].imass['Gypsum'] = self.inlet_solids # kg/hr. Mass of gypsum leaving in the solids outlet stream
        self.outs[0].imass['Nd'] = self.ins[0].imass['Nd'] # kg/hr. Mass of REEs leaving in the solids outlet stream
        self.outs[0].imass['U'] = self.ins[0].imass['U'] # kg/hr. Mass of radionuclides leaving in the solids outlet stream

        self.outs[1].imass['Water'] = self.ins[0].imass['Water'] # kg/hr. Amount of water leaving in the liquid outlet
        self.outs[1].imass['H2SO4'] = self.ins[0].imass['H2SO4'] # kg/hr. Amount of sulfuric acid leaving in the liquid outlet
        if self.Ln_desired == True: # includes REEs in the outlet stream if designated
            self.outs[1].imass['Nd'] = self.ins[0].imass['Nd'] # kg/hr. Amount of REE leavingin the liquid outlet

    _units = {
        'Total Filter Area': 'm2',
        'Number of Units': '#',
        'Filter Area per Unit': 'm2',
        'Total Air Suction': 'm3/min', 
        'Air Suction per Pump': 'm3/min',
        'Number of Pumps': '#'
    }

    def _design(self):
        D = self.design_results
        m2_to_f2 = 10.764 # conversion. 1 m2 = 10.764 ft2

        # filter
        self.max_area = 800/m2_to_f2 # m2. The largest size a vacuum rotary drum filter can be to use the cost correlations from Seider
        if self.filter_area < self.max_area: # if the size factor is smaller than the max size factor, then only cost one unit
            self.num_units = 1
            self.area = self.filter_area
        else: # if the size factor is larger than the max size factor, then split the filter area into several units and cost them
            self.num_units = self.filter_area/self.max_area
            self.area = self.max_area

        # vacuum pump
        air_suction = self.air_flow*self.area # m3/min. Air flow through vacuum system. 1 is middle of range from Hendriksson, Brandt - 2000 - focus on separation in the mining industry
        ft3min_to_m3min = 0.0283 # conversion. 1 ft3/min = 0.0283 m3/min
        max_air_suction = 350*ft3min_to_m3min # m3/min. The maximum flow at suction for a liquid-ring pump from Seider and Seader pg 595
        if air_suction > max_air_suction: # if the required suction is higher than the max suction from Seider, use multiple pumps
            num_pumps = air_suction/max_air_suction # number of pumps required
            suction = max_air_suction # m3/min. Suction of the individual pumps
        else: # if the required suction is less than the max from Seider, only use one pump
            num_pumps = 1 # number of pumps required
            suction = air_suction # m3/min. Suction of the individual pumps
        self.cp_pump = 8250*(suction/ft3min_to_m3min)**0.35*num_pumps # $. Purchase cost of the vacuum pumps. Seider and Seader pg 595

        D['Total Filter Area'] = self.filter_area
        D['Filter Area per Unit'] = self.area
        D['Number of Units'] = self.num_units
        D['Total Air Suction'] = self.filter_area
        D['Air Suction per Pump'] = suction
        D['Number of Pumps'] = num_pumps
        D['Purchase Cost of Vacuum Pumps'] = self.cp_pump
        

    def _cost(self):
        CEPCI = 1.6652 # CEPCI adjustment (2022/2006)
        m2_to_f2 = 10.764 # conversion. 1 m2 = 10.764 ft2
        Fm = 3.2 # material factor. from Seider and Seader pg 576. 
        self.baseline_purchase_costs['Rotary Drum Filter'] = \
            (np.exp(11.67-0.1905*np.log(self.area*m2_to_f2)+0.0554*np.log(self.area*m2_to_f2)**2)*self.num_units*Fm + self.cp_pump)*CEPCI # $. Purchase price for rotary drum filters from Seider pg 594
        
        self.power_utility.consumption = self.vacuum_energy*self.filter_area # kW based on range for power consumption of vacuum filters in Hendriksson, Brandt - 2000 - focus on separation in the mining industry



# %%
