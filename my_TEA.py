"""
Phosphogypsum Rare Earth Element Recovery System (PREERS)

Last Modified: 02/08/23

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""
import biosteam as bst
import qsdsan as qs
import numpy as np
import bst_TEA_modified

class my_TEA(bst_TEA_modified.TEA):
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, WC_over_FCI, 
                 maintenance, labor, membrane_cost
                 ): 
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor,
                         # Assume construction can be done within 1 year
                         construction_schedule=(1,),
                         # Assume no startup period
                         startup_months=0, startup_FOCfrac=0,
                         startup_VOCfrac=0, startup_salesfrac=0,
                         # Assume no financing
                         finance_interest=0, finance_years=0, finance_fraction=0,
                         # Assume no working capital
                         WC_over_FCI=0)
        
        # Create Attributes
        self.construction_schedule = construction_schedule
        self.WC_over_FCI = WC_over_FCI
        self.maintenance = maintenance
        self.labor = labor
        self.membrane_cost = membrane_cost
    
    # The abstract _DPI method should take installed equipment cost
    # and return the direct permanent investment. 
    def _DPI(self, installed_equipment_cost):
        # Lang factor already applied to installed_equipment_cost in bst_TEA_modified.py
        return installed_equipment_cost

    # The abstract _TDC method should take direct permanent investment
    # and return the total depreciable capital. 
    def _TDC(self, DPI):
        return (self.DPI ) # + contingencies_contractor_fees

    # The abstract _FCI method should take total depreciable capital
    # and return the fixed capital investment. 
    def _FCI(self, TDC):
        return (self.TDC + self.membrane_cost) # DPI, TDC, and FCI are all the same since the lang factor already includes the other costs of the capital investment # + land + royalties + startup
        # membrane cost included here because the cost of the unit is astronomically high. Using a factor to estimate installation cost would greatly overestimate the capital investment

    # The abstract _FOC method should take fixed capital investment
    # and return the fixed operating cost.
    def _FOC(self, TDC):
        # Calculated according to Seider - 4th edition product and process design principles pg 500
        # Labor Costs included in my_TEA.__init__ from calculation in systems.py

        # Maintenance Costs
        MWandB = self.maintenance*self.TDC # maintenance wages and benefits. For solid/liquid process (4.5% TDC)
        M_salaries_benefits = 0.25*MWandB # maintenance salaries and benefits
        M_supplies = 1*MWandB # maintenance supplies, materials, services and other
        M_overhead = 0.05*MWandB # maintenance overhead
        self.FOC_maintenance = MWandB + M_salaries_benefits + M_supplies + M_overhead # Total maintenance costs

        # Operating Overhead
        MandO_SWandB = self.labor + MWandB + M_salaries_benefits # maintenance and operations salary, wages, and benefits
        general_plant = 0.071*MandO_SWandB
        mechanical_dep_services = 0.024*MandO_SWandB
        employee_relations_dep = 0.059*MandO_SWandB
        business_services = 0.074*MandO_SWandB
        self.FOC_overhead = general_plant + mechanical_dep_services + employee_relations_dep + business_services

        # Other
        self.FOC_property_taxes_insurance = 0.02*self.TDC # property taxes and insurance

        # General Expenses
        selling_expense = 0.03*self.system.sales
        # direct_research = 0.048*self.system.sales
        # allocated_research = 0.005*self.system.sales
        # admin_expense = 0.02*self.system.sales
        # management_incentive = 0.0125*self.system.sales
        self.FOC_general_expenses = selling_expense # + direct_research + allocated_research + admin_expense + management_incentive not included due to this process being assumed NOAK and a waste treatment process

        return (self.labor + self.FOC_maintenance + self.FOC_overhead + self.FOC_property_taxes_insurance + self.FOC_general_expenses)