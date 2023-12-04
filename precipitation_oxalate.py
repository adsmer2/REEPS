import qsdsan as qs
import numpy as np

# Create the class
class precipitation_oxalate(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    under_solids_conc = 0.05, # weight fraction of underflow that is solids. Perry's 1999 Table 18-7 has examples of reasonable values
    res_time = 120, # mins. Residence time
    OA_uncertainty = 1 # multiplies feed rate of OA into the system for uncertainty/sensitivity
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes
        self.under_solids_conc = under_solids_conc
        self.res_time = res_time
        self.OA_uncertainty = OA_uncertainty
        
    _N_ins = 2
    _N_outs = 3

    def _run(self):
        # Calculate mass of oxalic acid needed for precipitation of all Lns
        MW_Ln = 144.24 # g/mol. Molecular weight of neodymium 
        MW_OA = 90.03 # g/mol. Molecular weight of oxalic acid
        MW_gypsum = 136.14 # # g/mol. Molecular weight of oxalic acid
        ratio_OA_Ln = 1.5 # the mol ratio of OA to Ln to form Ln2(OA)3
        ratio_OA_gypsum = 1 # the mol ratio of OA to gypsum to form Ca-OA complex
        solubility_gypsum = 0.2/100 # kg CaSO4/kg water
        self.OA_Ca = solubility_gypsum*self.ins[0].imass['water']/MW_gypsum*MW_OA*ratio_OA_gypsum *self.OA_uncertainty # kg/hr. stoichiometric oxalic acid required for coordinating just the Ca2+. Caclulated from the solubility of CaSO4 in solution
        self.OA_Ln = self.ins[0].imass['Nd']/MW_Ln*MW_OA*ratio_OA_Ln *self.OA_uncertainty # kg/hr stoichiometric oxalic acid required for coordinating just the Ln3+. The 1.5 is the mol ratio of OA to Ln
        excessOA_factor = 1.25 # an overconcentration factor used by industry to ensure high Ln recovery from Qi, Dezhi - 2018 - hydrometallurgy of REEs.
        self.OA_needed = (self.OA_Ca + self.OA_Ln)*excessOA_factor # (kg/hr) Qi, Dezhi - 2018 - hydrometallurgy of REEs. 1.25 is an overconcentration factor to increase Ln recovery. Notes in OneNote. 

        # Set Inlet Oxalic Acid
        self.ins[1].imass['OA'] = self.OA_needed # kg/hr. Oxalic acid feed

        # Define Outlet/Underflow Flows
        self.outs[0].imass['Nd'] = self.ins[0].imass['Nd'] # kg/hr. All REEs are precipitated
        complexation_adjust_factor = 88.02/90.03 # This factor is required to close the mass balance. It uses the MW of OA w/o and w/ Hydrogens, respectively, to account for the release of 2H atoms from OA during complexation to the REE atom. 
        self.outs[0].imass['OA'] = self.OA_Ln*complexation_adjust_factor # kg/hr. Mass of OA in the REE-oxalate precipitate
        
        solids = self.outs[0].imass['Nd'] + self.outs[0].imass['OA'] # kg/hr. Total mass of the REE-oxalate precipitate
        liquid = (1-self.under_solids_conc)*(solids/(self.under_solids_conc)) # kg/hr. Total mass of the liquid part of the underflow

        frac_sulf_in = self.ins[0].imass['H2SO4']/self.ins[0].imass['Water'] # wt fraction. fraction of liquid that comes from sulfuric acid
    
        # for simplicity, consider only OA in the solid the underflow
        self.outs[0].imass['Water'] = liquid*(1-frac_sulf_in) # kg/hr. mass of water in the underflow with the solids
        self.outs[0].imass['H2SO4'] = liquid*frac_sulf_in # kg/hr. mass of sulfuric acid in the underflow with the solids

        # Define Overflow Flows
        # Assumes 100% of solids goes to underflow
        self.outs[1].imass['Water'] = self.ins[0].imass['Water'] - self.outs[0].imass['Water'] # kg/hr. Mass of water in the overflow
        self.outs[1].imass['H2SO4'] = self.ins[0].imass['H2SO4'] - self.outs[0].imass['H2SO4'] # kg/hr. Mass of nitric acid in the overflow
        self.outs[1].imass['OA'] = self.OA_needed - self.OA_Ln # kg/hr the remaining soluble OA in wastewater assuming no precipitaiton of Ca-oxalate
        self.outs[1].imass['U'] = self.ins[0].imass['U'] # kg/hr of radionuclides in the overflow

        # account for hydrogen atoms that dissociate and become water soluble after complexation
        self.outs[2].imass['H+'] = self.OA_Ln - self.outs[0].imass['OA']  # mass of OA needed for stoich complexation - OA w/o hydrogens in the complex. The *90.03/88.019 the MW of OA w/ and w/o Hydrogens. Isolating the mass of hydrogen from the flow

    _units = {
        'Volume': 'm3',
        'Settling Area': 'm2',
        'Diameter': 'm',
        'Height': 'm',
        'OA for Ln':'kg/hr'
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
        D['OA for Calcium'] = self.OA_Ca*1.25
        D['OA for Ln'] = self.OA_Ln*1.25
        D['OA feed'] = self.OA_needed
        

    def _cost(self):
        m2_to_f2 = 10.7639 # 1 m2 = 10.7639 ft2
        Fm = 3.2 # material factor of 3.2
        CEPCI = 1.6652 # CEPCI adjustment to 2022 dollars 
        self.baseline_purchase_costs['OA Precipitation'] = \
            3360*(self.design_results['Settling Area']*m2_to_f2)**0.58*Fm*CEPCI # Purchase price with material factor of 3.2 from Seider and Seader pg 576 and CEPCI adjustment of 1.6652
        
        # Scale assuming the electricity usage is proportional to the volumetric flow rate
        self.power_utility.consumption = 1.2*(self.design_results['Diameter']/60) # kJ/hr. Power consumption for agitation. Assumed a proportional relationship to the Perry's thickener power consumption (12 kW for a thickener w/ dia = 60m from Perry's chemical engineers' handbook)


# %%
