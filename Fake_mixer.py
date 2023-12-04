import qsdsan as qs
import numpy as np

# Create the class
class Fake_mixer(qs.SanUnit):
    def __init__(self, ID= '', ins=None, outs=(), thermo=None, init_with='SanStream'
    ):
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
    _N_ins = 2
    _ins_size_is_fixed = False
    _N_outs = 1

    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins)

    _units = {
    }

    def _design(self):
        pass

    def _cost(self):
        pass