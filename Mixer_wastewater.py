import qsdsan as qs
import numpy as np

# Create the class
class Mixer_wastewater(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    # operating_price = 0.5 # Price for removal of organics is $0.5/kg inflation adjusted from Product and Process Design Principles, Synthesis, Analysis and Design, Third Edition - Seider et al
    operating_price = 0.0623 # $/m3 Operating cost for WWT inflation adjusted from Analysis, synthesis, and design of chemical processes, fifth edition - Richard Turton et al.

    ):  
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes
        self.operating_price = operating_price
        
        
    _N_ins = 3
    _N_outs = 1

    def _run(self):
        
        # Define outlet stream flowrates
        self.outs[0].imass['H2SO4'] = self.ins[0].imass['H2SO4'] + self.ins[1].imass['H2SO4'] + self.ins[2].imass['H2SO4']
        self.outs[0].imass['H2O'] = self.ins[0].imass['H2O'] + self.ins[1].imass['H2O'] + self.ins[2].imass['H2O'] 
        self.outs[0].imass['OA'] = self.ins[0].imass['OA'] + self.ins[1].imass['OA'] + self.ins[2].imass['OA']
        self.outs[0].imass['HNO3'] = self.ins[0].imass['HNO3'] + self.ins[1].imass['HNO3'] + self.ins[2].imass['HNO3']
        self.outs[0].imass['NaOH'] = self.ins[0].imass['NaOH'] + self.ins[1].imass['NaOH'] + self.ins[2].imass['NaOH']
        self.outs[0].imass['Na3PO4'] = self.ins[0].imass['Na3PO4'] + self.ins[1].imass['Na3PO4'] + self.ins[2].imass['Na3PO4']

    _units = {
        'OA removed':'kg/hr',
        'Wastewater to Treat':'gal/min'
    }

    def _design(self):
        D = self.design_results
        D['OA removed'] = self.outs[0].imass['OA']
        D['Wastewater to Treat'] = self.outs[0].F_vol*264.172/60
        

    def _cost(self):
        self.baseline_purchase_costs['Wastewater Treatment'] = 88000*(self.outs[0].F_vol*264.172/60)**0.64*1.6652/2 # purchase cost of WWTP assuming bare module factor of 2. Size factor in gal/min

        # self.add_OPEX = {'Wastewater Treatment': self.operating_price*self.outs[0].imass['OA']} # $/hr. 
        self.add_OPEX = {'Wastewater Treatment': self.operating_price*self.outs[0].F_vol} # $/hr
