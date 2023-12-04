import qsdsan as qs
import numpy as np

# Create the class
class PGmixer(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream',
    feedPG = 1000000, # kg/hr PG being fed to the process
    REEinPG = 0.5/100, # wt fraction of REE in the feed PG
    UinPG = 0.003159/100 # wt fraction of U in PG. From Liang et al – 2017 – Rare earths recovery and gypsum upgrade from Florida phosphogypsum
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes 
        self.feedPG = feedPG 
        self.REEinPG = REEinPG # fraction of REE in the feed PG
        self.UinPG = UinPG
        
    _N_ins = 1
    _N_outs = 1

    def _run(self):
        # Define inlet flowrates
        REEtotal = self.feedPG*self.REEinPG # kg/hr. Mass of REEs within the PG
        self.ins[0].imass['Nd'] = REEtotal # kg/hr. Set inlet flow rate of REEs within the PG feed

        Utotal = self.feedPG*self.UinPG # kg/hr. The total mass of radionuclides in the PG.
        self.ins[0].imass['U'] = Utotal # kg/hr. Set inlet flow rate of radionuclides within the PG feed

        gypsum = self.feedPG - REEtotal - Utotal # kg/hr. the mass of gypsum in the PG feed stream based on the composition data from Liang et al - 2017
        self.ins[0].imass['Gypsum'] = gypsum # kg/hr. Set the inlet flow of gypsum within the PG feed

        # Define outlet stream flowrates
        self.outs[0].imass['Gypsum'] = self.ins[0].imass['Gypsum']
        self.outs[0].imass['Nd'] = self.ins[0].imass['Nd']
        self.outs[0].imass['U'] = self.ins[0].imass['U']

    _units = {
    }

    def _design(self):
        pass

    def _cost(self):
        pass