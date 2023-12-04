import qsdsan as qs
import numpy as np
import math

# Create the class
class resuspension_tank(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    dil_ratio = 4, # ratio of liquid to solids after resuspension
    
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes 
        self.dil_ratio = dil_ratio

        self.F_BM = {'resuspension tank': 2} # need to confirm this

    _N_ins = 3
    _N_outs = 1

    def _run(self):
        # Calculate mass of acid needed to resuspend solids
        mol_Ln2O3 = self.ins[0].imol['Nd2O3'] # kmol/hr Ln2O3 in feed
        # Ln2O3 + 6HNO3 -> 2Ln(NO3)3 + 3H2O. It's 6 mol HNO3 per mol Ln2O3 for solution to remain neutral.
        mass_acid = mol_Ln2O3*6*63.01284/0.68 # kg/hr of 68% (w/w) HNO3. Nitric acid most commercially available at 68% (w/w) purity. MW 63.01
        mass_water = mass_acid*(1-0.68) # kg/hr mass of water in the purchased acid
        mass_water_dissolution = mol_Ln2O3*3*15.9994 # kg/hr. from rxn eqn, 1 mol Ln2O3 yields 3 mol H2O. Use 15.999 instead of 18.02 since hydrogen's are counted with the nitric acid.
        mass_excess_water = mass_acid*self.dil_ratio - mass_acid # kg/hr total mass of liquid in the system

        # Define required DI water for inlet stream
        # self.ins[1].imass['HNO3'] = mass_acid # kg/hr
        self.ins[1].imol['HNO3'] = mol_Ln2O3*6/0.68
        self.ins[2].imass['H2O'] = mass_excess_water # kg/hr
        
        # Define outlet stream flowrates
        self.outs[0].imol['Nd'] = mol_Ln2O3*2 # kmol/hr. from rxn eqn, 1 mol Ln2O3 yields 2 mol Ln(NO3)3
        # self.outs[0].imass['HNO3'] = mass_acid - mass_water # kg/hr. Due to the dilution of the HNO3 feed, only pure HNO3 is in this stream for accurate volume calcs. Mass balance was confirmed to close
        self.outs[0].imol['HNO3'] = mol_Ln2O3*6
        self.outs[0].imass['H2O'] = mass_water + mass_excess_water + mass_water_dissolution # kg/hr

    _units = {
        'Volume': 'm3',
        'Diameter': 'm',
        'Thickness': 'm',
        'Weight SS': 'kg', 
        'Storage Residence Time': 'days',
        'Vol. Flow of Lixiviant In': 'm3/hr',
        'Volume of Tank': 'gal',
        'Cost of Storage Tank': '$'
    }

    def _design(self):
        D = self.design_results
        HRT = 0.25 # hydraulic residence time in hours
        self.outlet_flowrate = qs.SanStream('outlet_flowrate')
        self.outlet_flowrate.mix_from(self.outs)
        w_vol = (self.outlet_flowrate.F_vol)*HRT # working volume (m3)
        tot_vol = w_vol / 0.8 # total volume assuming 80% working volume (m3)
        dia = (16*tot_vol/3.14)**(1/3) # diameter with dia/height=4 (m)
        H = dia/4
        D['Volume'] = tot_vol
        D['Diameter'] = dia

        # Calculate wall thickness based on pressure vessel formuala in Seider, Seader
        # Pd = 10 + 14.7 # design pressure psi
        # stress = 18750 # for 304 stainless steel (psi) https://www.nickelinstitute.org/media/1842/types304and304lstainlesssteelsforlowtemperatureservice_328_.pdf
        # weld = 0.9 # fractional weld efficiency 0.85-1 from https://www.pveng.com/joint-efficiency/
        # self.tp = Pd*dia/(stress*weld - 1.2*Pd)*0.0254 # pressure vessel wall thickness (m) Seider and Seader pg 575
        self.tp = 0.009525 # minimum allowed thickness (m) 
        row_steel = 8000 # kg/m3
        self.steel_m = (3.14*H*((self.tp + dia)**2 - dia**2)/4 + 3.14*self.tp*dia**2/4)*row_steel # mass of steel for annulus plus tank bottom. Must be 1;000 < W < 920;000 lb:
        D['Thickness'] = self.tp
        D['Weight SS'] = self.steel_m

        # Nitric Acid Storage Tank
        # ------------------------
        # assume process water for lixiviant will be piped in
        # Use a cone roof tank 
        days_storage = 7 # days of sulfuric acid stored assuming weekly deliveries. Varys from 1 week to 1 month. from Seider and Seader pg 588
        Qlix_in = self.ins[1].F_vol # m3/hr fresh sulfuric acid feedstock flowrates
        Vstorage = Qlix_in*264.172*24 *days_storage # gal of stored sulfuric acid
        self.cost_storage = 265*(Vstorage**0.51) *1.6652*1.7 # $. Material factor 1.7 SS-304 for nitric acid resistance. Adjusted by CEPCI 1.6652 (2022/2006) from Seider and Seader pg 595 

        D['--Storage Tank--'] = ''
        D['Storage Residence Time'] = days_storage
        D['Vol. Flow of Lixiviant In'] = Qlix_in
        D['Volume of Tank'] = Vstorage
        D['Cost of Storage Tank'] = self.cost_storage
        

    def _cost(self):
        self.baseline_purchase_costs['resuspension tank'] = \
            math.exp(8.9552-0.23330*math.log(self.steel_m*2.20462) + 0.04333*math.log(self.steel_m*2.20462)**2)*1.7*1.6652 + self.cost_storage # Purchase price with material factor of 1.7 for SS-304 from Seider and Seader pg 576 
        
        # Mixing tank needs impellar utility consumption
        # self.power_utility.consumption = 1.2*(self.design_results['Diameter']/60) 