import qsdsan as qs
import numpy as np

# Create the class
class precipitation_NaPhosphate(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    under_solids_conc = 0.05, # weight fraction of underflow that is solids. Perry's 1999 Table 18-7 has examples of reasonable values
    res_time = 120, # mins. residence time
    NaOH_uncertainty = 1, # multiplies feed rate of NaOH into the system for uncertainty/sensitivity
    Na3PO4_uncertainty = 1, # multiplies feed rate of Na3PO4 into the system for uncertainty/sensitivity
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes
        self.under_solids_conc = under_solids_conc
        self.res_time = res_time
        self.NaOH_uncertainty = NaOH_uncertainty
        self.Na3PO4_uncertainty = Na3PO4_uncertainty
        
    _N_ins = 4
    _N_outs = 2

    def _run(self):
        # Calculate mass of NaOH and Na3PO4 needed for precipitation of all U
        # Visual Minteq simulation gave the pH (~9) and concentration of precipitant required for full precipitation (>99.99%) of U as Na-Autunite (Na2(UO2)2(PO4)2)
        ratio_sodphos_U = 25 # the mol ratio of Na3PO4 to U to form (Na2(UO2)2(PO4)2)
        MW_sodphos = 163.94 # g/mol
        MW_NaOH = 39.997 # g/mol
        MW_U = 238.026 # g/mol
        MW_H2SO4 = 98.079 # g/mol

        self.sodphos_needed = self.ins[0].imass['U']/MW_U*ratio_sodphos_U  *MW_sodphos *self.Na3PO4_uncertainty # kg/hr sodium phosphate. keeping the ratio of 25 mol phosphate to mol U as in the VisualMinteq simulation
        # Neutralize solution
        ratioNaOH_neutralize = 2 # ratio of NaOH needed to neutralize sulfuric acid according to H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O 
        sodphos_neutralize = 1 # amount of hydrogen that sodium phosphate neutralizes at pH 9 according to H+ + Na2PO4,1- -> HNa2PO4 @ pH 9
        self.NaOH_needed = (self.ins[0].imass['H2SO4']/MW_H2SO4*ratioNaOH_neutralize *MW_NaOH - self.sodphos_needed/MW_sodphos*sodphos_neutralize *MW_NaOH)*self.NaOH_uncertainty # kg/hr. Total NaOH needed to get to pH 9

        self.sodphos_U = self.sodphos_needed/ratio_sodphos_U # kg/hr. Amount of sodium phosphate actually in the precipitate Na2(UO2)2(PO4)2

        # Set Inlet sodium phosphate and NaOH
        self.ins[1].imass['Na3PO4'] = self.sodphos_needed # kg/hr. trisodium phosphate feed
        self.ins[2].imass['NaOH'] = self.NaOH_needed # kg/hr. sodium hydroxide feed

        # Define Outlet/Underflow Flows
        self.outs[0].imass['U'] = self.ins[0].imass['U'] # kg/hr. Mass of radionuclides leaving in the underflow solids. Visual minteq predicted >99.99% (assume 100%) precipitation
        self.outs[0].imass['Na3PO4'] = self.sodphos_U # kg/hr. Trisodium phosphate leaving in the underflow solids
        solids = self.outs[0].imass['U'] + self.outs[0].imass['Na3PO4'] # kg/hr. Total solids flow rate in underflow. 
        liquid = (1-self.under_solids_conc)*(solids/(self.under_solids_conc)) # kg/hr. Total liquid flow rate based on the specified underflow solids content

        frac_sulf_in = (self.ins[0].imass['H2SO4'] + self.ins[3].imass['H2SO4'])/(self.ins[0].imass['Water'] + self.ins[3].imass['Water']) # wt fraction. fraction of liquid that is sulfuric acid

        self.outs[0].imass['Water'] = liquid*(1-frac_sulf_in) # kg/hr. mass of water in the underflow with the solids
        self.outs[0].imass['H2SO4'] = liquid*frac_sulf_in # kg/hr. mass of sulfuric acid in the underflow with the solids

        # Define Outlet/Overflow Flows
        # Assumes 100% of solids goes to underflow
        self.outs[1].imass['Water'] = self.ins[0].imass['Water'] + self.ins[3].imass['Water'] - self.outs[0].imass['Water'] # kg/hr. Mass of water in the overflow
        self.outs[1].imass['H2SO4'] = self.ins[0].imass['H2SO4'] + self.ins[3].imass['H2SO4'] - self.outs[0].imass['H2SO4'] # kg/hr. Mass of sulfuric acid in the overflow
        self.outs[1].imass['OA'] = self.ins[0].imass['OA'] # kg/hr. Mass of oxalic acid in the overflow. Assume all OA goes to wastewater.
        self.outs[1].imass['NaOH'] = self.ins[2].imass['NaOH'] # kg/hr. Mass of sodium hydroxide in the overflow
        self.outs[1].imass['Na3PO4'] = self.sodphos_needed - self.sodphos_U # kg/hr. Mass of the remaining Na3PO4 in wastewater

    _units = {
        'Volume': 'm3',
        'Settling Area': 'm2',
        'Diameter': 'm',
        'Height': 'm',
        'Na3PO4 feed':'kg/hr'
    }

    def _design(self):
        D = self.design_results
        HRT = self.res_time/60 # hrs. hydraulic residence time based on leaching time
        self.outlet_flowrate = qs.SanStream('outlet_flowrate') # Initialize a new stream
        self.outlet_flowrate.mix_from(self.outs) # Make this new stream a mixture of all outlet flows
        w_vol = (self.outlet_flowrate.F_vol)*HRT # m3. Working volume assuming
        tot_vol = w_vol / 0.8 # m3. Total volume of each tank assuming 80% working volume
        dia = (24*tot_vol/3.14)**(1/3) # m. Diameter of each vessel
        D['Volume'] = tot_vol
        D['Settling Area'] = 3.14*(dia/2)**2
        D['Diameter'] = dia
        D['Height'] = H = dia/6
        D['Na3PO4 feed'] = self.sodphos_needed
        D['NaOH feed'] = self.NaOH_needed
        

    def _cost(self):
        m2_to_f2 = 10.7639 # 1 m2 = 10.7639 ft2
        Fm = 3.2 # material factor of 3.2
        CEPCI = 1.6652 # CEPCI adjustment to 2022 dollars 
        self.baseline_purchase_costs['OA Precipitation'] = \
            3360*(self.design_results['Settling Area']*m2_to_f2)**0.58*Fm*CEPCI # Purchase price from Seider and Seader pg 576 
        
        # Scale assuming the electricity usage is proportional to the volumetric flow rate
        self.power_utility.consumption = 12*(self.design_results['Diameter']/60) # kJ/hr. Power consumption for agitation. Assumed a proportional relationship to the Perry's thickener power consumption (12 kW for a thickener w/ dia = 60m from Perry's chemical engineers' handbook)


# %%
