"""
Phosphogypsum Rare Earth Element Recovery System (PREERS)

Last Modified: 02/07/23

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""
import qsdsan as qs
import numpy as np

# Create the leaching class
class fired_heater_oxalate_pure(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    temp = 850, # temperature needed to decompose oxalates for all Lns (degrees C) from Qi, Dezhi - 2018
    time = 1.5, # time in furnace (hr) heat for 1.5-2 hrs from Qi, Dezhi - 2018
    # heat for 1.5-2 hrs from Qi, Dezhi - 2018
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes 
        self.temp = temp
        self.time = time
        
    _N_ins = 1
    _N_outs = 3

    def _run(self):
        self.inlet_flowrate = qs.SanStream('inlet_flowrate')
        self.inlet_flowrate.mix_from(self.ins)

        # Calculate volume of furnace and air flow rate
        s_vol = (self.inlet_flowrate.F_vol)*self.time # volume of solids in furnace (m3)
        self.ratio_air = 1 # ratio of air to solids
        self.tot_vol = s_vol*(1 + self.ratio_air) # include extra volume for air
        cp_air = 1.006 # kJ/kg/K @ 20oC https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html
        self.v_air_flow = self.inlet_flowrate.F_vol*self.ratio_air # m3/hr
        self.density_air = 1.204 # kg/m3 @ 20oC https://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html

        # Heat Duty Calculation
        self.heat_eff = 0.8 # efficiency of fired heater

        # Cp of inlet stream ~1.097 kJ/kg/K
        self.solids_HD = self.inlet_flowrate.Cp/1000*1000*self.inlet_flowrate.F_mass*(self.temp-20) # heat duty for solids in kJ/hr
        self.air_HD = self.v_air_flow*self.density_air*cp_air*(self.temp-20) # heat duty for air in kJ/hr

        # oxalate combustion reaction    Ln2(C2O4)3 · 10H2O(l) → Ln2O3 + 3CO + 3CO2 + 10H2O(g) consider the water as an additional heating cost
        mol_Ln = self.ins[0].imol['Nd']*1000 # mol/hr of lanthanide in the feed
        mol_Ln2OA3 = mol_Ln/2 # 2 mol Ln per mol of the precipitated complex Ln2OA3
        mol_water = 10*mol_Ln2OA3 # mol/hr of water in feed
        mol_CO2 = 3*mol_Ln2OA3 # mol/hr of CO2 produced
        mol_CO = 3*mol_Ln2OA3 # mol/hr of CO produced
        mol_Nd2O3 = mol_Ln2OA3 # mol/hr of Ln calcined to its oxide
        cp_water = 75.38 # J/mol/K https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Type=JANAFL&Plot=on
        dHvap_water = 40.65 # kJ/mol
        self.water_HD = mol_water*cp_water*(self.temp-20)/1000 + dHvap_water*mol_water # kJ/hr
        dHc = -242.9 # kJ/mol heat of combustion of oxalic acid from https://webbook.nist.gov/cgi/cbook.cgi?ID=C144627&Mask=2
        self.OA_comb_HD = mol_Ln*3/2 *dHc # kJ/hr

        self.heat_duty = (self.solids_HD + self.air_HD + self.water_HD + self.OA_comb_HD)/self.heat_eff # heat duty in kJ/hr
        
        # Natural Gas and Air Required
        LHV_nat_gas = 47100 # kJ/kg Seider and Seader and https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html
        MW_nat_gas = 19
        mol_nat_gas = self.heat_duty/LHV_nat_gas/MW_nat_gas # kmol/hr
        #using complete combustion of methane, get required oxygen -> air
        O_content_air = 0.21 # content of oxygen in the air
        self.m_air = mol_nat_gas*2*32/O_content_air # kg/hr
        self.m_CO2_nat_gas = mol_nat_gas*44.01 # kg/hr CO2 produced by nat gas combustion

        # Solid REO Product
        self.outs[0].imol['Nd2O3'] = mol_Nd2O3/1000 # kmol/hr

        # Combustion emissions
        self.outs[1].imass['CO2'] = mol_CO2*44.01/1000 # kg/hr CO2 from heating with NG accounted for in utlities
        self.outs[2].imass['CO'] = mol_CO*28.01/1000 # kg/hr
        # self.outs[3].imass['Water'] = mol_water*18.02/1000 # kg/hr. no impact on process so disregarded (water content of solids wasn't modeled so it doesn't enter the mass balance)

    _units = {
        'Heat Duty Total': 'kJ/hr @ 80 pcnt efficiency',
        'Volume': 'm3',
        'air/solids volume ratio': '-',
        'Air Flow Rate': 'kg/hr',
        'Heat Duty Solids': 'kJ/hr',
        'Heat Duty Air': 'kJ/hr',
        'Heat Duty Water': 'kJ/hr',
        'Efficiency': '%',
        'Air for Combustion': 'kg/hr',
        'Heat Duty OA Combustion': 'kJ/hr'
    }

    def _design(self):
        D = self.design_results
              
        D['Volume'] = self.tot_vol
        D['air/solids volume ratio'] = self.ratio_air
        D['Air Flow Rate'] = self.v_air_flow*self.density_air

    
        D['Efficiency'] = self.heat_eff*100
        D['Heat Duty Total'] = self.heat_duty
        D['Heat Duty Solids'] = self.solids_HD
        D['Heat Duty Air'] = self.air_HD
        D['Heat Duty Water'] = self.water_HD
        D['Heat Duty OA Combustion'] = self.OA_comb_HD

        D['Air for Combustion'] = self.m_air
        D['CO2 by combustion of nat gas'] = self.m_CO2_nat_gas
        D['CO2 by combustion of oxalate'] = self.outs[1].imass['CO2']
        

    def _cost(self):
        self.baseline_purchase_costs['Fired Heater'] = np.exp(0.32325 + 0.766*np.log(self.design_results['Heat Duty Total']*0.948))*1.6652 # converted kJ to btu and CEPCI adjusted. Cost eqn from Seider design book

        # Heat Utility Cost
        self.add_OPEX = {'Natural Gas 3': 7.534*10**(-6)*self.heat_duty} # 7.534E-6 $/kJ * heat duty kJ/hr = $/hr. NG price for industry July 2022 from EIA