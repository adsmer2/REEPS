"""
Phosphogypsum Rare Earth Element Recovery System (PREERS)

Last Modified: 10/17/22

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""

# Import packages
from logging import raiseExceptions
from chaospy import distributions 
import qsdsan as qs
import numpy as np

#
def create_model_leaching(sys, fununit):
    model = qs.Model(sys)
    param = model.parameter
    metric = model.metric

    tea = sys.TEA
    lca = sys.LCA

    flowsheet = qs.Flowsheet.flowsheet.default
    fs_stream = flowsheet.stream
    fs_unit = flowsheet.unit
    
    
    # Leaching Unit
    U1 = fs_unit.U1

    baseline = U1.time
    dist = distributions.Uniform(lower=15, upper=205)
    @param(name='Leaching Time', element='U1', kind='coupled', units='mins', baseline=baseline, distribution=dist)
    def set_leachingTime(i):
        U1.time = i

    baseline = U1.temp
    dist = distributions.Uniform(lower=21, upper=69)
    @param(name='Leaching Temperature', element='U1', kind='coupled', units='deg C', baseline=baseline, distribution=dist)
    def set_leachingTemp(i):
        U1.temp = i

    baseline = U1.acidConc
    dist = distributions.Uniform(lower=0.5, upper=9.5)
    @param(name='Acid Concentration', element='U1', kind='coupled', units='wt pcnt acid', baseline=baseline, distribution=dist)
    def set_leaching_acidConc(i):
        U1.acidConc = i

    baseline = U1.solventRatio
    dist = distributions.Uniform(lower=2.1, upper=6.9)
    @param(name='Solvent/Solid Ratio', element='U1', kind='coupled', units='Ratio of Liquid/Solid', baseline=baseline, distribution=dist)
    def set_leaching_solventRatio(i):
        U1.solventRatio = i

    # ----------------------------------------------------------------------------------------
    # Metrics
    # ----------------------------------------------------------------------------------------

    # Economic
    # -----------
    @metric(name='NPV15', units='MM USD', element='TEA')
    def get_NPV():
        return tea.NPV/1000000

    @metric(name='IRR', units='%', element='TEA')
    def get_IRR():
        return tea.solve_IRR()
        
    @metric(name='MSP', units='USD/kg REO', element='TEA')
    def get_MSP():
        return tea.solve_price(fs_stream.Ln2O3)

    # # Environmental
    # # -----------
    # @metric(name='Global Warming', units=f'kg CO2-Eq/kg {fununit}', element='LCA')
    # def get_annual_GWP100():
    #     if fununit == 'REO':
    #         return lca.total_impacts['GWP100']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['GWP100']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['GWP100']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Urban Land Occupation', units=f'square meter-year/kg {fununit}', element='LCA')
    # def get_annual_ULOP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['ULOP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['ULOP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['ULOP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Agricultural Land Occupation', units=f'square meter-year/kg {fununit}', element='LCA')
    # def get_annual_ALOP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['ALOP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['ALOP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['ALOP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Freshwater Ecotoxicity', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    # def get_annual_FETPinf():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['FETPinf']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['FETPinf']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['FETPinf']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Human Toxicity', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    # def get_annual_HTPinf():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['HTPinf']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['HTPinf']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['HTPinf']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Metal Depletion', units=f'kg Fe-Eq/kg {fununit}', element='LCA')
    # def get_annual_MDP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['MDP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['MDP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['MDP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Fossil Depletion', units=f'kg oil-Eq/kg {fununit}', element='LCA')
    # def get_annual_FDP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['FDP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['FDP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['FDP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Freshwater Eutrophication', units=f'kg P-Eq/kg {fununit}', element='LCA')
    # def get_annual_FEP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['FEP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['FEP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['FEP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Terrestrial Acidification', units=f'kg SO2-Eq/kg {fununit}', element='LCA')
    # def get_annual_TAP100():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['TAP100']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['TAP100']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['TAP100']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Particulate Matter Formation', units=f'kg PM10-Eq/kg {fununit}', element='LCA')
    # def get_annual_PMFP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['PMFP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['PMFP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['PMFP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Ionising Radiation', units=f'kg U235-Eq/kg {fununit}', element='LCA')
    # def get_annual_IRP_HE():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['IRP_HE']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['IRP_HE']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['IRP_HE']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Marine Ecotoxicity', units=f'kg 1,4-DB-Eq/kg {fununit}', element='LCA')
    # def get_annual_METPinf():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['METPinf']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['METPinf']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['METPinf']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Terrestrial Ecotoxicity', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    # def get_annual_TETPinf():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['TETPinf']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['TETPinf']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['TETPinf']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Natural Land Transformation', units=f'square meter/kg {fununit}', element='LCA')
    # def get_annual_NLTP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['NLTP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['NLTP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['NLTP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Photochemical Oxidant Formation', units=f'kg NMVOC-Eq/kg {fununit}', element='LCA')
    # def get_annual_POFP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['POFP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['POFP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['POFP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Marine Eutrophication', units=f'kg N-Eq/kg {fununit}', element='LCA')
    # def get_annual_MEP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['MEP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['MEP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['MEP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Ozone Depletion', units=f'kg CFC-11-Eq/kg {fununit}', element='LCA')
    # def get_annual_ODPinf():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['ODPinf']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['ODPinf']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['ODPinf']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    # @metric(name='Water Depletion', units=f'm3 water-Eq/kg {fununit}', element='LCA')
    # def get_annual_WDP():
    #     if fununit == 'Ln':
    #         return lca.total_impacts['WDP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
    #     elif fununit == 'PG':
    #         return lca.total_impacts['WDP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
    #     elif fununit == 'none':
    #         return lca.total_impacts['WDP']
    #     else:
    #         raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')

    return model