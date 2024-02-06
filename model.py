"""
Rare Earth Element recovery from Phosphogypsum System (REEPS)

The Pennsylvania State University
Chemical Engineering Department
S2D2 Lab (Dr. Rui Shi)
@author: Adam Smerigan
"""

# Import packages
from logging import raiseExceptions
from chaospy import distributions 
import qsdsan as qs
import numpy as np

__all__ = ('create_model',)

#
def create_model(sys, fununit, parameter, target):
    model = qs.Model(sys)
    param = model.parameter
    metric = model.metric

    tea = sys.TEA
    lca = sys.LCA

    flowsheet = qs.Flowsheet.flowsheet.default
    fs_stream = flowsheet.stream
    fs_unit = flowsheet.unit

    # ----------------------------------------------------------------------------------------
    # Technological Parameters
    # ----------------------------------------------------------------------------------------
    def set_technological_params():
        # Leaching
        # -------------
        U1 = fs_unit.U1

        baseline = U1.time # 200
        dist = distributions.Triangle(lower=190, midpoint=baseline, upper=205)
        @param(name='Leaching Time (U1)', element='U1', kind='coupled', units='mins', baseline=baseline, distribution=dist)
        def set_leachingTime(i):
            U1.time = i

        baseline = U1.temp # 47
        dist = distributions.Triangle(lower=baseline-0.05*baseline, midpoint=baseline, upper=baseline+0.05*baseline) # assume 5% potential error from experiment
        @param(name='Leaching Temperature (U1)', element='U1', kind='coupled', units='deg C', baseline=baseline, distribution=dist)
        def set_leachingTemp(i):
            U1.temp = i

        baseline = U1.acidConc # 3.2
        dist = distributions.Triangle(lower=baseline-0.05*baseline, midpoint=baseline, upper=baseline+0.05*baseline)
        @param(name='Acid Concentration (U1)', element='U1', kind='coupled', units='wt pcnt acid', baseline=baseline, distribution=dist)
        def set_leaching_acidConc(i):
            U1.acidConc = i

        baseline = U1.solventRatio # 2.75
        dist = distributions.Triangle(lower=baseline-0.05*baseline, midpoint=baseline, upper=baseline+0.05*baseline)
        @param(name='Solvent to Solid Ratio (U1)', element='U1', kind='coupled', units='Ratio of Liquid/Solid', baseline=baseline, distribution=dist)
        def set_leaching_solventRatio(i):
            U1.solventRatio = i

        baseline = U1.LOverV
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='overflow to underflow ratio (U1)', element='U1', kind='coupled', units='mass fraction', baseline=baseline, distribution=dist)
        def set_leaching_LOverV(i):
            U1.LOverV = i

        # Precipitation 
        # -------------
        # P1 Precipitation Oxalate
        baseline = fs_unit.P1.OA_uncertainty
        dist = distributions.Triangle(lower=baseline*0.75, midpoint=baseline, upper=baseline*1.25)
        @param(name='Oxalate Feed (P1)', element='P1', kind='coupled', units='kg/hr', baseline=baseline, distribution=dist)
        def set_P1_oxalate_feed(i):
            fs_unit.P1.OA_uncertainty = i

        baseline = fs_unit.P1.res_time
        dist = distributions.Uniform(lower=60, upper=180)
        @param(name='Residence Time (P1)', element='P1', kind='coupled', units='min', baseline=baseline, distribution=dist)
        def set_P1_res_time(i):
            fs_unit.P1.res_time = i

        # P2 Precipitation Oxalate Pure
        baseline = fs_unit.P2.OA_uncertainty
        dist = distributions.Triangle(lower=baseline*0.85, midpoint=baseline, upper=baseline*1.15)
        @param(name='Oxalate Feed (P2)', element='P2', kind='coupled', units='kg/hr', baseline=baseline, distribution=dist)
        def set_P2_oxalate_feed(i):
            fs_unit.P2.OA_uncertainty = i

        baseline = fs_unit.P2.res_time
        dist = distributions.Uniform(lower=60, upper=180)
        @param(name='Residence Time (P2)', element='P2', kind='coupled', units='min', baseline=baseline, distribution=dist)
        def set_P2_res_time(i):
            fs_unit.P2.res_time = i

        # P3 Precipitation Na3PO4
        baseline = fs_unit.P3.NaOH_uncertainty
        dist = distributions.Triangle(lower=baseline*0.75, midpoint=baseline, upper=baseline*1.25)
        @param(name='Sodium Hydroxide Feed (P3)', element='P3', kind='coupled', units='kg/hr', baseline=baseline, distribution=dist)
        def set_P3_NaOH_uncertainty(i):
            fs_unit.P3.NaOH_uncertainty = i

        baseline = fs_unit.P3.Na3PO4_uncertainty
        dist = distributions.Triangle(lower=baseline*0.75, midpoint=baseline, upper=baseline*1.25)
        @param(name='Na3PO4 Feed (P3)', element='P3', kind='coupled', units='kg/hr', baseline=baseline, distribution=dist)
        def set_P3_Na3PO4_uncertainty(i):
            fs_unit.P3.Na3PO4_uncertainty = i

        baseline = fs_unit.P3.res_time
        dist = distributions.Uniform(lower=60, upper=180)
        @param(name='Residence Time (P3)', element='P3', kind='coupled', units='min', baseline=baseline, distribution=dist)
        def set_P3_res_time(i):
            fs_unit.P3.res_time = i

        # Filters
        # -------------
        # Filter_l
        baseline = fs_unit.F1.vacuum_energy
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Vacuum Energy (F1)', element='F1', kind='coupled', units='kW/m2', baseline=baseline, distribution=dist)
        def set_F1_vacuum_energy(i):
            fs_unit.F1.vacuum_energy = i

        baseline = fs_unit.F2.vacuum_energy
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Vacuum Energy (F2)', element='F2', kind='coupled', units='kW/m2', baseline=baseline, distribution=dist)
        def set_F2_vacuum_energy(i):
            fs_unit.F2.vacuum_energy = i

        baseline = fs_unit.F1.air_flow
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Air Flow (F1)', element='F1', kind='coupled', units='m3/m2 filter area/min', baseline=baseline, distribution=dist)
        def set_F1_air_flow(i):
            fs_unit.F1.air_flow = i

        baseline = fs_unit.F2.air_flow
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Air Flow (F2)', element='F2', kind='coupled', units='m3/m2 filter area/min', baseline=baseline, distribution=dist)
        def set_F2_air_flow(i):
            fs_unit.F2.air_flow = i

        # Black box
        # -------------
        S1 = fs_unit.S1

        baseline = fs_unit.S1.immobilization_density
        dist = distributions.Triangle(lower=baseline*0.1, midpoint=baseline, upper=baseline*10)
        @param(name='Immobilization Density (S1)', element='S1', kind='coupled', units='mmol/L adsorbent bed', baseline=baseline, distribution=dist)
        def set_S1_capacity(i):
            fs_unit.S1.capacity = i

        baseline = fs_unit.S1.cycle_time
        dist = distributions.Triangle(lower=2, midpoint=baseline, upper=24)
        @param(name='Cycle Time (S1)', element='S1', kind='coupled', units='hrs', baseline=baseline, distribution=dist)
        def set_S1_cycle_time(i):
            fs_unit.S1.cycle_time = i

        baseline = fs_unit.S1.membrane_lifetime
        dist = distributions.Triangle(lower=1, midpoint=baseline, upper=20)
        @param(name='Membrane Lifetime (S1)', element='S1', kind='coupled', units='years', baseline=baseline, distribution=dist)
        def set_S1_cycle_time(i):
            fs_unit.S1.membrane_lifetime = i

        # baseline = fs_unit.S1.pressure
        # dist = distributions.Triangle(lower=baseline*0.1, midpoint=baseline, upper=baseline*10)
        # @param(name='Pressure Drop (S1)', element='S1', kind='coupled', units='Pa', baseline=baseline, distribution=dist)
        # def set_S1_cycle_time(i):
        #     fs_unit.S1.pressure = i

        if target == 'no':
            baseline = S1.recovery
            dist = distributions.Uniform(lower=0.9, upper=1)
            @param(name='REE Recovery (S1)', element='S1', kind='coupled', units='fraction', baseline=baseline, distribution=dist)
            def set_S1_recovery(i):
                S1.recovery = i

            baseline = fs_unit.S1.capacity
            dist = distributions.Triangle(lower=0.0005, midpoint=baseline, upper=5)
            @param(name='Adsorbent Capacity (S1)', element='S1', kind='coupled', units='mol/L adsorbent', baseline=baseline, distribution=dist)
            def set_S1_capacity(i):
                fs_unit.S1.capacity = i

        elif target == 'REE Recovery (S1)':
            baseline = S1.recovery
            dist = distributions.Uniform(lower=0.2, upper=1)
            @param(name='REE Recovery (S1)', element='S1', kind='coupled', units='fraction', baseline=baseline, distribution=dist)
            def set_S1_recovery(i):
                S1.recovery = i

        elif target == 'Adsorbent Capacity (S1)':
            baseline = fs_unit.S1.capacity
            dist = distributions.Uniform(lower=0.002, upper=0.01)
            @param(name='Adsorbent Capacity (S1)', element='S1', kind='coupled', units='mol/L adsorbent', baseline=baseline, distribution=dist)
            def set_S1_capacity(i):
                fs_unit.S1.capacity = i

        else:
            raise RuntimeError('in function "create model" argument "target" must be either "no" or a prespecified parameter')
    
    # ----------------------------------------------------------------------------------------
    # Contextual Parameters
    # ----------------------------------------------------------------------------------------
    def set_contextual_params():
        # TEA Parameters
        # --------------------
        # Interest Rate
        baseline = tea.IRR
        dist = distributions.Uniform(lower=0.1, upper=0.2)
        @param(name='Interest Rate', element='TEA', kind='isolated', units='-', baseline=baseline, distribution=dist)
        def set_interest_rate(i):
            tea.IRR = i

        # Lang Factor
        baseline = tea.lang_factor
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Lang Factor', element='TEA', kind='isolated', units='-', baseline=baseline, distribution=dist)
        def set_lang_factor(i):
            tea.lang_factor = i

        # Operating Days
        baseline = tea.operating_days
        dist = distributions.Triangle(lower=365*0.8, midpoint=baseline, upper=365*0.95)
        @param(name='Operating Days', element='TEA', kind='isolated', units='days', baseline=baseline, distribution=dist)
        def set_operating_days(i):
            tea.operating_days = i

        # Income Tax
        baseline = tea.income_tax
        dist = distributions.Uniform(lower=0.21, upper=0.325)
        @param(name='Income Tax Rate', element='TEA', kind='isolated', units='-', baseline=baseline, distribution=dist)
        def set_income_tax(i):
            tea.income_tax = i

        # Labor Cost
        baseline = tea.labor
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Labor', element='TEA', kind='isolated', units='USD/year', baseline=baseline, distribution=dist)
        def set_labor(i):
            tea.labor = i
        
        # Stream Prices
        # ----------------------
        # Ln2O3 product price
        baseline = fs_stream.Ln2O3.price
        dist = distributions.Triangle(lower=baseline*0.7, midpoint=baseline, upper=baseline*1.3)
        @param(name='REO Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_Ln2O3_price(i):
            fs_stream.Ln2O3.price = i

        # Gypsum product price
        baseline = fs_stream.gypsum.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Gypsum Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_gypsum_price(i):
            fs_stream.gypsum.price = i

        # H2SO4 price
        baseline = fs_stream.lixiviant_acid.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='H2SO4 Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_H2SO4_Price(i):
            fs_stream.lixiviant_acid.price = i

        # OA price
        baseline = fs_stream.OA_feed.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Oxalic Acid Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_OA_feed_Price(i):
            fs_stream.OA_feed.price = i
            fs_stream.OA_feed_2.price = i

        # NaOH price
        baseline = fs_stream.NaOH_feed.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='NaOH Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_NaOH_feed_Price(i):
            fs_stream.NaOH_feed.price = i

        # Water price
        baseline = fs_stream.lixiviant_water.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Process Water Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_lixiviant_water_Price(i):
            fs_stream.lixiviant_water.price = i
            fs_stream.RS_water.price = i

        # Na3PO4 price
        baseline = fs_stream.Na3PO4_feed.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='Na3PO4 Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_Na3PO4_feed_Price(i):
            fs_stream.Na3PO4_feed.price = i

        # HNO3 price
        baseline = fs_stream.HNO3_feed.price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2)
        @param(name='HNO3 Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_HNO3_feed_Price(i):
            fs_stream.HNO3_feed.price = i

        # WWT operating cost
        baseline = fs_unit.WT.operating_price
        dist = distributions.Triangle(lower=baseline*0.5, midpoint=baseline, upper=baseline*1.5)
        @param(name='Wastewater Treatment Price', element='TEA', kind='isolated', units='$/kg', baseline=baseline, distribution=dist)
        def set_WT_operating_price(i):
            fs_unit.WT.operating_price = i

        # Electricity Price
        baseline = qs.PowerUtility.price
        dist = distributions.Triangle(lower=.06, midpoint=baseline, upper=0.16) # https://www.eia.gov/electricity/data/browser/#/topic/7?agg=1,0&geo=vvg&endsec=2&freq=M&start=200101&end=202305&ctype=columnchart&ltype=pin&rtype=s&maptype=0&rse=0&pin=
        @param(name='Electricity Price', element='TEA', kind='isolated', units='$/kWh', baseline=baseline, distribution=dist)
        def set_PowerUtility_price(i):
            qs.PowerUtility.price = i

        # Peptide Price (S1)
        baseline = fs_unit.S1.peptide_price
        dist = distributions.Triangle(lower=baseline*0.3, midpoint=baseline, upper=baseline*1.2) 
        @param(name='Peptide Price', element='TEA', kind='isolated', units='$/g peptide', baseline=baseline, distribution=dist)
        def set_S1_peptide_price(i):
            fs_unit.S1.peptide_price = i

        # Resin Price (S1)
        baseline = fs_unit.S1.resin_price
        dist = distributions.Triangle(lower=baseline*0.8, midpoint=baseline, upper=baseline*1.2) 
        @param(name='Resin Price', element='TEA', kind='isolated', units='$/L resin', baseline=baseline, distribution=dist)
        def set_S1_resin_price(i):
            fs_unit.S1.resin_price = i

        # Characterization Factor
        # -----------------------
        # Heating from Natural Gas
        # baseline = sys.heatNG_item.CFs['GWP']
        # dist = distributions.LogNormal()
        # @param(name='Stainless steel GWP', element='LCA', kind='isolated', units='kg CO2/kg', baseline=baseline, distribution=dist)
        # def set_heatNG_GWP(i):
        #     sys.heatNG_item.CFs['GWP'] = i

    if parameter == 'all':
        set_technological_params()
        set_contextual_params()
    elif parameter == 'technological':
        set_technological_params()
    elif parameter == 'contextual':
        set_contextual_params()
    else:
        raise RuntimeError(f'In create_model(sys, fununit, parameter), parameter={parameter} is not "technological" or "contextual". Please define as one of these two.')


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
        return tea.solve_IRR()*100
        
    @metric(name='MSP', units='USD/kg REO', element='TEA')
    def get_MSP():
        return tea.solve_price(fs_stream.Ln2O3)

    # Environmental
    # -----------
    @metric(name='Acidification Terrestrial', units=f'kg SO\u2082-Eq/kg {fununit}', element='LCA')
    def get_annual_TAP():
        if fununit == 'Ln':
            return lca.total_impacts['TAP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['TAP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['TAP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Climate Change', units=f'kg CO\u2082-Eq/kg {fununit}', element='LCA')
    def get_annual_GWP100():
        if fununit == 'Ln':
            return lca.total_impacts['GWP1000']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['GWP1000']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['GWP1000']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Ecotoxicity Freshwater', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    def get_annual_FETP():
        if fununit == 'Ln':
            return lca.total_impacts['FETP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['FETP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['FETP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Ecotoxicity Marine', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    def get_annual_METP():
        if fununit == 'Ln':
            return lca.total_impacts['METP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['METP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['METP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Ecotoxicity Terrestrial', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    def get_annual_TETP():
        if fununit == 'Ln':
            return lca.total_impacts['TETP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['TETP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['TETP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Energy Resources', units=f'kg oil-Eq/kg {fununit}', element='LCA')
    def get_annual_FFP():
        if fununit == 'Ln':
            return lca.total_impacts['FFP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['FFP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['FFP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Eutroph. Freshwater', units=f'kg P-Eq/kg {fununit}', element='LCA')
    def get_annual_FEP():
        if fununit == 'Ln':
            return lca.total_impacts['FEP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['FEP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['FEP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Eutroph. Marine', units=f'kg P-Eq/kg {fununit}', element='LCA')
    def get_annual_MEP():
        if fununit == 'Ln':
            return lca.total_impacts['MEP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['MEP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['MEP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Human Toxicity Carc.', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    def get_annual_HTPc():
        if fununit == 'Ln':
            return lca.total_impacts['HTPc']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['HTPc']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['HTPc']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Human Toxicity N-carc.', units=f'kg 1,4-DCB-Eq/kg {fununit}', element='LCA')
    def get_annual_HTPnc():
        if fununit == 'Ln':
            return lca.total_impacts['HTPnc']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['HTPnc']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['HTPnc']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Ionising Radiation', units=f'kBq Co-60-Eq/kg {fununit}', element='LCA')
    def get_annual_IRP():
        if fununit == 'Ln':
            return lca.total_impacts['IRP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['IRP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['IRP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Land Use', units=f'm\u00B2*a crop-Eq/kg {fununit}', element='LCA')
    def get_annual_LOP():
        if fununit == 'Ln':
            return lca.total_impacts['LOP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['LOP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['LOP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Meterial Resources', units=f'kg Cu-Eq/kg {fununit}', element='LCA')
    def get_annual_SOP():
        if fununit == 'Ln':
            return lca.total_impacts['SOP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['SOP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['SOP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Ozone Depletion', units=f'kg CFC-11-Eq/kg {fununit}', element='LCA')
    def get_annual_ODPinfinite():
        if fununit == 'Ln':
            return lca.total_impacts['ODPinfinite']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['ODPinfinite']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['ODPinfinite']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Particulate Matter', units=f'kg PM2.5-Eq/kg {fununit}', element='LCA')
    def get_annual_PMFP():
        if fununit == 'Ln':
            return lca.total_impacts['PMFP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['PMFP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['PMFP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Photochemical Ox. Human Health', units=f'kg NOx-Eq/kg {fununit}', element='LCA')
    def get_annual_HOFP():
        if fununit == 'Ln':
            return lca.total_impacts['HOFP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['HOFP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['HOFP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Photochemical Ox. Ecosystems', units=f'kg NOx-Eq/kg {fununit}', element='LCA')
    def get_annual_EOFP():
        if fununit == 'Ln':
            return lca.total_impacts['EOFP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['EOFP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['EOFP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')
        
    @metric(name='Water Use', units=f'cubic meter/kg {fununit}', element='LCA')
    def get_annual_WCP():
        if fununit == 'Ln':
            return lca.total_impacts['WCP']/(lca.lifetime*fs_stream.Ln2O3.F_mass*tea.operating_days*24)
        elif fununit == 'PG':
            return lca.total_impacts['WCP']/(lca.lifetime*fs_stream.rawPG.F_mass*tea.operating_days*24)
        elif fununit == 'none':
            return lca.total_impacts['WCP']
        else:
            raise NameError(f'For model metric, {fununit} is not "Ln" or "PG"')

    return model