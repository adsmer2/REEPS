import qsdsan as qs
import numpy as np

# Create the class
# Assume the membrane adsorbption system creates a saleable pure product
class black_box(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    recovery = 0.99, # recovery of lanthanides (% as decimal)
    capacity = 0.00577, # mol/L agarose adsorbent beads. From ACS Cent. Sci. 2021, 7, 1798−1808
    # capacity = 0.025, # mol/L agarose adsorbent beads. Target Value
    # recovery = 1, # recovery of lanthanides (% as decimal). Target Value
    immobilization_density = 2.47, # mmol/L adsorbent bed. From ACS Cent. Sci. 2021, 7, 1798−1808
    cycle_time = 4, # hrs. The time it takes for one column to stop, start, and desorb REEs.
    peptide_price = 0.5, # $/g biomolecule. from multiple sources of protein and peptide data
    resin_price = 45, # $/L adsorbent resin. From https://samcotech.com/how-much-does-it-cost-to-buy-maintain-and-dispose-of-ion-exchange-resins/
    membrane_lifetime = 10, # years. How long before the membrane must be replaced?
    pressure = 101325*5 # Pa. Pressure drop through membrane system. Taken as the high end of ultrafiltration in Richard D. Noble (1987) An Overview of Membrane Separations, Separation Science and Technology, 22:2-3, 731-743
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes 
        self.recovery = recovery
        self.capacity = capacity
        self.immobilization_density = immobilization_density
        self.cycle_time = cycle_time
        self.peptide_price = peptide_price
        self.resin_price = resin_price
        self.membrane_lifetime = membrane_lifetime
        self.pressure = pressure
        
    _N_ins = 1
    _N_outs = 2

    def _run(self):
        # This unit is assumed to make ~99.99% or greater purity individual REE streams
        # Regeneration isn't considered here since experimental data is lacking for adsorption systems for this application
        # Define outlet stream flowrates
        self.outs[0].imass['Nd'] = self.ins[0].imass['Nd']*self.recovery
        self.outs[0].imass['H2O'] = self.ins[0].imass['H2O']
        self.outs[0].imass['HNO3'] = self.ins[0].imass['HNO3']

        # Have a separate outlet for REEs that are not recovered so that the mass balance closes
        self.outs[1].imass['Nd'] = self.ins[0].imass['Nd']*(1-self.recovery)

    _units = {
        'Flow Rate': 'm3/hr',
        'Ln Flow In': 'mol/hr',
        'Ln Recovery': '%',
        'Membrane Lifetime':'years',
        'Membrane Capital Cost': '$',
        'Membrane Operating Cost': '$/year',
        'Peptide Needed': 'g'
    }

    def _design(self):
        D = self.design_results
        
        # Estimate Capital/Operating Cost
        # -----------------------
        # Calculate adsorbent cost based on ion exchange resin data
        Ln_in = self.ins[0].imol['Nd']*1000 # mol/hr
        resin_needed = Ln_in/self.capacity*self.cycle_time *2 # L adsorbent resin. Need two adsorption columns to allow for regeneration of the other column
        MW_peptide = 11800 # g/mol molecular weight of the LanM protein
        peptide_needed = resin_needed*self.immobilization_density/1000*MW_peptide # g peptide
        resin_cost = resin_needed*self.resin_price # $
        peptide_cost = peptide_needed*self.peptide_price # $ 

        self.membrane_cost = resin_cost + peptide_cost # cost of the first adsorbent purchase (capital cost)

        plant_lifetime = 30 # years
        replacements = (plant_lifetime - self.membrane_lifetime)/self.membrane_lifetime # years of plant life remaining after first adsorbent needs replacement/adsorbent lifetime. The number of times new adsorbent needs to be bought after the initial purchase
        self.annual_membrane_cost = self.membrane_cost*replacements/plant_lifetime # operating cost of adsorbent replacement annualized over the plant lifetime
        
        D['Flow Rate'] = self.ins[0].F_vol
        D['Ln Flow In'] = self.ins[0].imol['Nd']*1000
        D['Ln Recovery'] = self.recovery*100
        D['Peptide Needed'] = peptide_needed
        D['Peptide Cost'] = peptide_cost 
        D['Resin Cost'] = resin_cost
        D['Membrane Lifetime'] = self.membrane_lifetime
        D['Membrane Capital Cost'] = self.membrane_cost
        D['Membrane Operating Cost'] = self.annual_membrane_cost
        

    def _cost(self):
        self.baseline_purchase_costs['black box'] = 0 # The capital cost here is made zero to avoid using the Lang factor for this piece of equipment. The membrane_cost is added when calling my_TEA in systems.py (adding directly to TDC)
    
        # OPEX (membrane replacements)
        self.add_OPEX = {'Peptide Membrane Replacements': self.annual_membrane_cost/8760} # cost of replacing adsorbent ($/hr)
