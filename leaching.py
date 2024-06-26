"""
Rare Earth Element recovery from Phosphogypsum System (REEPS)

The Pennsylvania State University
Chemical Engineering Department
S2D2 Lab (Dr. Rui Shi)
@author: Adam Smerigan
"""
# %%
import numpy as np
import qsdsan as qs
from scipy import interpolate

# Create the leaching class
class LeachingSulfuric(qs.SanUnit):
    def __init__(self, ID='U1', ins=None, outs=(), thermo=None, init_with='SanStream',
    acidConc = 2.8, # wt %. acid concentration
    temp = 47, # degrees C. Leaching temperature 
    time = 200, # min. Leaching time 
    solventRatio = 2.75, # solvent to solid mass flow ratio in the leaching unit
    LOverV = 0.2, # underflow to overflow mass flow ratio
    feedPG = 1000000, # kg/hr. The mass flow rate of PG being fed to the system (also the plant capacity)
    REEinPG = 0.005 # wt fraction. REE content in the feed PG
    ):
        # Some Assumed Process Parameters
        self.Y_np1 = 0 # amount of REE remaining in recycled solvent feed (mass REE/mass solvent)
        self.X_N = 0.00001 # amount of REE leaving in underflow with solids (mass REE/mass solvent)

        # Some standard codes you need to include for all subclasses of `SanUnit`
        qs.SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        # These are the unique attributes of `LeachingSulfuric`
        self.acidConc = acidConc
        self.temp = temp
        self.time = time
        self.solventRatio = solventRatio
        self.LOverV = LOverV
        self.feedPG = feedPG
        self.REEinPG = REEinPG
        
    _N_ins = 4
    _N_outs = 2

    def _run(self):
        PG, lixiviant_water, lixiviant_acid, rec_lixiviant = self.ins # specify the names of inlet streams
        leachate, underflow = self.outs # specify the names of outlet streams

        # Add in Experimental Data to calculate leaching efficiency (LE)
        # data from Liang et al – 2017 – Rare earths recovery and gypsum upgrade from Florida phosphogypsum
        LEbase = 0.43 # base case leaching efficiency
        
        LEbase_temp = 0.43 # leaching efficiency for Liang's reported optimal temperature
        xtemp = [20, 30, 40, 50, 60, 70] # x data (temperature oC)
        ytemp = [0.33, 0.4, 0.42, 0.43, 0.39, 0.34] # leaching efficiency (fraction of PG REEs that are leached/total PG REEs)
        tck_temp = interpolate.splrep(xtemp, ytemp, k=1, s=0) # straight line interpolation of the experimental data so that we can find LE at any value of xdata

        LEbase_conc = 0.4 # leaching efficiency for Liang's reported optimal acid concentration
        xconc = [0, 1, 2.5, 5, 7.5, 10] # x data (acid concentration in wt %)
        yconc = [0.01, 0.18, 0.36, 0.4, 0.39, 0.35] # leaching efficiency (fraction of PG REEs that are leached/total PG REEs)
        tck_conc = interpolate.splrep(xconc, yconc, k=1, s=0) # straight line interpolation of the experimental data so that we can find LE at any value of xdata

        LEbase_time = 0.43 # leaching efficiency for Liang's reported optimal leaching time
        xtime = [10,20,30,45,63,76,90,103,116,135,150,165,180,190,210]
        ytime = [0.19,0.3,0.35,0.38,0.4,0.41,0.42,0.425,0.43,0.4325,0.435,0.435,0.435,0.4375,0.44] # leaching efficiency (fraction of PG REEs that are leached/total PG REEs)
        tck_time = interpolate.splrep(xtime, ytime, k=1, s=0) # straight line interpolation of the experimental data so that we can find LE at any value of xdata

        LEbase_solventRatio = 0.42 # leaching efficiency for Liang's chosen optimal solvent to solid ratio
        xsolventRatio = [2, 3, 4, 5, 6, 7] # x data (solvent/solid ratio)
        ysolventRatio = [0.32, 0.38, 0.42, 0.46, 0.50, 0.53] # leaching efficiency (fraction of PG REEs that are leached/total PG REEs)
        tck_solventRatio = interpolate.splrep(xsolventRatio, ysolventRatio, k=1, s=0) # straight line interpolation of the experimental data so that we can find LE at any value of xdata

        # Create errors if the model is running outside the bounds of the given experimental data
        if not 20 <= self.temp <= 70:
            raise RuntimeError('`temp` must be within [20, 70], '
                                f'the provided value {self.temp} is outside this range.')
        if not 0 <= self.acidConc <= 10:
            raise RuntimeError('`acidConc` must be within [0, 10], '
                                f'the provided value {self.acidConc} is outside this range.')
        if not 10 <= self.time <= 210:
            raise RuntimeError('`time` must be within [10, 210], '
                                f'the provided value {self.time} is outside this range.')
        if not 2 <= self.solventRatio <= 7:
            raise RuntimeError('`solventSolidRatio` must be within [2, 7], '
                                f'the provided value {self.solventRatio} is outside this range.')

        # This section calculates the difference in LE from the base case optimal reported by Liang et al (0.43) due to changing values of temperature, acid conc, leaching time, and solvent/solid ratio. 
        # ----------------------------------------------------
        # Changes in LE that are 'positive' indicate that the LE is lower than the base case LE and will be subtracted from the base case LE when calculating the final LE.
        LEtemp = LEbase_temp - interpolate.splev(self.temp, tck_temp, der=0).tolist() # difference in LE from the base case based on changing the temperature in this model (self.temp)
        LEconc = LEbase_conc - interpolate.splev(self.acidConc, tck_conc, der=0).tolist() # difference in LE from the base case based on changing the acid concentration in this model (self.acidconc)
        LEtime = LEbase_time - interpolate.splev(self.time, tck_time, der=0).tolist() # difference in LE from the base case based on changing the leaching time in this model (self.time)
        LEsolventRatio = LEbase_solventRatio - interpolate.splev(self.solventRatio, tck_solventRatio, der=0).tolist() # difference in LE from the base case based on changing the Solvent/solid ratio in this model (self.solventRatio)
        self.leachingEFF = LEbase - (LEtemp + LEconc + LEtime + LEsolventRatio) # sum the changes in LE from the base case for each parameter and subtract it from the base case LE to get the final LE at the desired parameter values

        # Calculate the required lixiviant to leach a given amount of PG (plant capacity)
        lixNeeded = self.feedPG*self.solventRatio # kg/hr. Total mass of lixiviant being used
        fracSulfuric = self.acidConc/100 # fraction of lixiviant that is sulfuric acid
        lixiviant_acid.imass['H2SO4'] = fracSulfuric*lixNeeded - self.ins[3].imass['H2SO4'] # kg/hr. Total sulfuric acid entering the system. Considers how much is being recycled in the second term
        lixiviant_water.imass['H2O'] = (1-fracSulfuric)*lixNeeded - self.ins[3].imass['Water'] # kg/hr. Total water entering the system. Considers how much is being recycled in the second term
        
        # Set flowrate of PG stream in and the amount of REE to be leached
        REEtotal = self.feedPG*self.REEinPG # kg/hr. Mass of REEs within the PG
        REEleached = REEtotal*self.leachingEFF #kg/hr. Total mass flow of REE being leached out of the PG in the leaching stage

        # Complete the Mass Balance for ideal leaching stages
        # --------------------------------------------
        # Based off the design equations for a countercurrent continuous leaching system w/ ideal leaching stages in Seader, Henley, Roper - Separation Process Principles 3rd edition
        # See the assumptions in the SI of the paper
        V_1 = lixNeeded # kg/hr. Flow rate of overflow liquid flowing into leaching stage
        V_L = V_1 - self.LOverV*V_1 # kg/hr. Leachate flow rate out of leaching stage
        L_L = V_1 - V_L # kg/hr. Underflow liquid flow rate out of leaching stage
        Y_L = (REEleached - L_L*self.X_N)/V_L # wt fraction. REE content in leachate in overflow (exiting the leaching stage)
        Y_1 = (V_L*Y_L + L_L*Y_L - REEleached)/V_1 # wt fraction. REE content entering leaching stage in overflow

        # make these values callable in the design section to calculate the # of theoretical stages
        self.Y_L = Y_L # REE content in the overflow leaving the leaching stage (mass REE/mass solvent)
        self.Y_1 = Y_1 # REE content in the overflow entering the leaching stage (mass REE/mass solvent)

        # Set outlet flows 
        # -------------------
        leachate.imass['Nd'] = Y_L*V_L # kg/hr. REEs leaving the leaching unit in the leachate (overflow)
        leachate.imass['H2SO4'] = fracSulfuric*V_L # kg/hr. Sulfuric acid leaving the leaching unit in the leachate (overflow)
        leachate.imass['H2O'] = (1-fracSulfuric)*V_L # kg/hr. Process water leaving teh leaching unit in the leachate (overflow)
        leachate.imass['Gypsum'] = self.ins[0].imass['Gypsum']*0.3 # kg/hr. In thickener some solids are in the overflow unlike we assumed for ideal leaching stages in the mass balance. 30% of the solids are assumed to remain in the overflow through operational control. Assume no gypsum dissociates/associates in leaching
        leachate.imass['U'] = self.ins[0].imass['U']*0.7 # kg/hr. 70% of U is leached assuming fine particles of PG. From Liang et al – 2017 – Rare earths recovery and gypsum upgrade from Florida phosphogypsum

        underflow.imass['Nd'] = self.X_N*L_L + (REEtotal-REEleached) # kg/hr. REEs leaving the leaching unit in the underflow
        underflow.imass['H2SO4'] = fracSulfuric*L_L # kg/hr. Sulfuric acid leaving the leaching unit in the underflow 
        underflow.imass['H2O'] = (1-fracSulfuric)*L_L # kg/hr. Process water leaving teh leaching unit in the underflow 
        underflow.imass['Gypsum'] = self.ins[0].imass['Gypsum']*0.7 # kg/hr. Read leachate solids comment
        underflow.imass['U'] = self.ins[0].imass['U']*0.3 # kg/hr. See above comment for leachate
        
 

    _units = {
        'Washing Stages': '#',
        'Total Stages': '#',
        'Volume': 'm3',
        'Settling Area': 'm2',
        'Diameter': 'm',
        'Height': 'm',
        'Thickness': 'm',
        'Steel': 'kg',
        'Heat Duty': 'kJ/h',
        'Leaching Efficiency': '%',
        'Volumetric Flowrate': 'm3/hr',
        'Pump Power Consumption': 'kW',
        'Storage Residence Time': 'days',
        'Vol. Flow of Lixiviant In': 'm3/hr',
        'Volume of Tank': 'gal',
        'Cost of Storage Tank': '$'
    }

    def _design(self):
        D = self.design_results # make design results easier to call
        stage_eff = 0.9 # fraction. Account for inefficiencies on leaching stages. Assumption
        self.nStages = np.log10(((self.X_N - self.Y_np1)/(self.Y_L - self.Y_1)))/np.log10(self.LOverV)/stage_eff # number of actual washing stages. Calculate number of actual washing stages based on McCabe-Smith algebraic method
        HRT = self.time/60 # hr. hydraulic residence time based on leaching time
        self.outlet_flowrate_leach = qs.SanStream('outlet_flowrate_leach') # initiate a new stream
        self.outlet_flowrate_leach.mix_from(self.outs) # assign all outlets to this stream for use in calculating total volume in the unit
        w_vol = (self.outlet_flowrate_leach.F_vol)*HRT/(self.nStages + 1) # m3. Working volume of each leaching/washing stage (thickener)
        tot_vol = w_vol / 0.8 # m3. total volume assuming 80% working volume
        dia = (24*tot_vol/3.14)**(1/3) # m. Assume aspect ratio of D/H=6
        D['Leaching Efficiency'] = self.leachingEFF*100
        D['Washing Stages'] = self.nStages
        D['Total Stages'] = (self.nStages + 1)
        D['Volume'] = tot_vol # m3
        D['Settling Area'] = 3.14*(dia/2)**2 # m2
        D['Diameter'] = dia # m
        D['Height'] = H = dia/6 # aspect ratio of 6
        D['Volumetric Flowrate'] = self.outlet_flowrate_leach.F_vol # m3/hr

        # Heat Utility
        self.inlet_flowrate_leach = qs.SanStream('inlet_flowrate') # Initiate a new stream
        self.inlet_flowrate_leach.mix_from(self.ins) # Make this new stream a mixture of all inlet flows
        self.heat_duty = self.inlet_flowrate_leach.Cp/1000*1000*self.inlet_flowrate_leach.F_mass*(self.temp-20) # kJ/hr. Heat duty to heat the lixiviant and gypsum to the specified temperature
        D['Heat Duty'] = self.heat_duty # kJ/hr

        # Centrifugal Pump and Motor for underflow
        # ----------------------------------------------
        # Calculations from pg 561 and 562 in Seider and Seader
        m3_to_gal = 264.172 # conversion. 1 m3 = 264.172 gal
        Qtank = self.outlet_flowrate_leach.F_vol/(self.nStages+1)*m3_to_gal/60 # gal/min. Flowrate out of each tank
        m_to_ft = 3.281 # conversion. 1 m = 3.281 ft
        head = H*1.5 *m_to_ft # ft. pump head 
        rowS = 2300 # kg/m3. Density of gypsum solids
        rowL = 1000 # kg/m3. Density of water
        cS = 0.7 # weight fraction of solids in underflow
        cL = 0.3 # weight fraction of liquids in underflow
        lbsgal_to_kgm3 = 119.8 # conversion. 1 lb/gal = 119.8 kg/m3
        density = 1/(cS/rowS+cL/rowL) / lbsgal_to_kgm3 # lbs/gal. Density of slurry 
        spec_grav = density/rowL # specific gravity of slurry
        effp = -0.316 + 0.24015*np.log(Qtank) - 0.01199*(np.log(Qtank))**2 # fraction. pump efficiency
        pump_break = Qtank*head*spec_grav/3960/effp # hp. Pump break power
        effm = 0.8 + 0.0319*np.log(pump_break) - 0.00182*np.log(pump_break)**2 # fraction. motor efficiency
        hp_to_kW = 0.746 # conversion. 1 hp = 0.746 kW
        self.powerc = Qtank*head*density/33000/effp/effm*hp_to_kW # kW. Power consumption 

        Spump = Qtank*(head)**0.5 # scaling factor for a centrifugal pump
        if Qtank >= 3500:
            Ft_pump = 2 # cost factor for type of pump
        elif Qtank < 3500 and Qtank >= 900:
            Ft_pump = 1.5
        elif Qtank < 900 and Qtank > 0: 
            Ft_pump = 1
        else:
            raise RuntimeError('flow rate is outside of reasonable limits')
        Fm_pump = 3 # pump material factor. assuming 3 for special application with slurry
        self.Cpump = np.exp(9.7171 - 0.6019*np.log(Spump) + 0.0519*np.log(Spump)**2) * Ft_pump*Fm_pump *(self.nStages+1) # purchase cost of centrifugal pump. from Seider and Seader pg 576. Adjusted by CEPCI 1.6652 (2022/2006)

        Ft_motor = 1.35 # motor type 
        self.Cmotor = np.exp(5.8259 + 0.13141*np.log(self.powerc) + 0.053255*np.log(self.powerc)**2 + 0.028628*np.log(self.powerc)**3 - 0.0035549*np.log(self.powerc)**4) *Ft_motor *(self.nStages+1) # purchase cost of motor with enclosed type. from Seider and Seader pg 576. Adjusted by CEPCI 1.6652 (2022/2006)

        D['--Pump--'] = ''
        D['Qtank'] = Qtank
        D['Pb'] = pump_break
        D['effp'] = effp
        D['head'] = head
        D['Cost of Pump'] = self.Cpump
        D['Cost of Motor'] = self.Cmotor
        D['Pump Power Consumption'] = self.powerc

        # Belt Conveyor for feed PG
        # ---------------------
        elevate = H # m. Height the solids need to travel vertically
        L = elevate/np.sin(np.pi/6) # length of belt in m. Calculated assuming a 30 degree incline up the vertical rise of the leaching tank. Assume W = 1 middle of reasonable range from Seider. Assume H is in reasonable range by adjusting rate of conveying
        W = 1 # width of belt in m
        self.WL = W*39.37*L*3.28 # width in inches and length in ft for cost equation

        kghr_to_lbs = 0.000612395 # conversion. 1 kg/hr = that number lb/sec
        self.conv_P = 0.00058*((self.ins[0].F_mass*kghr_to_lbs)**0.82)*L*m_to_ft*hp_to_kW + 0.00182*self.ins[0].F_mass*kghr_to_lbs*elevate*m_to_ft*hp_to_kW # kW. Power consumption of conveyor. from Seider and Seader pg 588
        self.conv_C = 21.5*self.WL # $. Purchase cost of belt conveyor from from Seider and Seader pg 595

        D['--conveyor--'] = ''
        D['lenth of belt'] = L
        D['width of belt'] = W
        D['height to elevate solids'] = elevate
        D['Cost of Conveyor'] = self.conv_C
        D['Conveyor Power Consumption'] = self.conv_P

        # Storage tank cost for sulfuric acid
        # -------------------------
        # assume process water for lixiviant will be piped in
        # Use a cone roof tank 
        days_storage = 7 # days of sulfuric acid stored assuming weekly deliveries. Varys from 1 week to 1 month. from Seider and Seader pg 588
        Qlix_in = self.ins[2].F_vol # m3/hr fresh sulfuric acid feedstock flowrates
        Vstorage = Qlix_in*m3_to_gal*24 *days_storage # gal of stored sulfuric acid. Max volume 20 million gallons from Seider and Seader pg 589 
        Fm_storagetank = 3.2 # material factor for the storage tank for sulf acid resistance
        self.cost_storage = 265*(Vstorage**0.51) *Fm_storagetank # $. Purchase cost of a storage tank with a conical roof. from Seider and Seader pg 595 

        D['--Storage Tank--'] = ''
        D['Storage Residence Time'] = days_storage
        D['Vol. Flow of Lixiviant In'] = Qlix_in
        D['Volume of Tank'] = Vstorage
        D['Cost of Storage Tank'] = self.cost_storage

        

    def _cost(self):
        m2_to_f2 = 10.7639 # 1 m2 = 10.7639 ft2
        Fm_leaching = 3.2 # Material factor of 3.2 (carpenter 20CB-3) from Seider and Seader pg 576. 
        CEPCI = 1.6652 # CEPCI adjustment (2022/2006)
        self.baseline_purchase_costs['Thickener'] = \
            (3360*(self.design_results['Settling Area']*m2_to_f2)**0.58*(self.nStages+1)*Fm_leaching    + self.Cpump + self.Cmotor + self.conv_C + self.cost_storage)*CEPCI # $. Purchase price of the entire leaching unit system
        
        # Scale assuming the electricity usage is proportional to the volumetric flow rate
        self.power_utility.consumption = 12*(self.design_results['Diameter']/60) + self.powerc + self.conv_P # Power consumption of the unit. Agitation estimated from 12 kW for a thickener w/ dia = 60m from Perry's chemical engineers' handbook

        # Heat Utility Cost
        price_NG = 7.534*10**(-6) # $/kJ. Price of NG. NG price for industry July 2022 from EIA
        self.add_OPEX = {'Natural Gas': price_NG*self.heat_duty} # $/hr. Hourly cost of NG use for heating

        
        

        

