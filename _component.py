"""
Phosphogypsum Rare Earth Element Recovery System (PREERS)

Last Modified: 02/01/23

The Pennsylvania State University
Chemical Engineering Department
@author: Adam Smerigan
"""
# %%
import qsdsan as qs

__all__ = ('create_components', )


def create_components():
#  Setup inlet flows
    Calcium_sulfate_dihydrate = qs.Chemical(ID='Calcium_sulfate_dihydrate', 
        CAS=10101-41-4, 
        search_db=True, 
        phase='s', 
        rho=2320, # kg/m3      https://pubchem.ncbi.nlm.nih.gov/compound/24928
        Cp=186.2/172.171, # J/g      http://chemister.ru/Database/properties-en.php?dbid=1&id=552
        Hf=-2021.1*1000) # J/mol   http://chemister.ru/Database/properties-en.php?dbid=1&id=552

    Nd2O3 = qs.Chemical(ID='Nd2O3', 
        # CAS=1313-97-9, 
        search_db=False,
        phase='s',
        MW = 336.48, # https://pubchem.ncbi.nlm.nih.gov/compound/Neodymium-oxide
        rho=7240, # kg/m3 http://www.chemspider.com/Chemical-Structure.140158.html       
        Cp=111.3/172.171, # J/g https://en.wikipedia.org/wiki/Neodymium(III)_oxide          
        Hf=-1807.9*1000) # J/molhttps://pubs.acs.org/doi/10.1021/je60039a030
    
    Na3PO4 = qs.Chemical(ID='Na3PO4', 
        CAS=7601-54-9, 
        search_db=True,
        phase='s',
        rho=1620, # kg/m3 https://www.chemspider.com/Chemical-Structure.22665.html
        Cp=665/163.939, # J/g http://chemister.ru/Database/properties-en.php?dbid=1&id=781         
        Hf=-5480*1000) # J/mol http://chemister.ru/Database/properties-en.php?dbid=1&id=781
    
    Hion = qs.Chemical(ID='H+', 
        CAS=12385-13-6, 
        search_db=True,
        phase='l',
        rho= 89.9, # kg/m3
        Cp=14.3, # J/g/K https://www.engineeringtoolbox.com/hydrogen-d_976.html
        Tb=373.12, # K Assumed the same for water since this is soluble H+ in water    
        # Hf=1530 # J/mol https://webbook.nist.gov/cgi/inchi?ID=C12385136&Mask=20#Ion-Energetics
        ) 

    # cmps_default = qs.Components.load_default()
    H2O_chem = qs.Chemical('H2O')
    H2O = qs.Component.from_chemical(ID='H2O', chemical=H2O_chem, particle_size='Soluble', degradability='Undegradable', organic=False)
    Hion = qs.Component.from_chemical(ID='H+', chemical=Hion, particle_size='Soluble', degradability='Undegradable', organic=False)
    Neodymium = qs.Component('Nd', search_ID='7440-00-8', particle_size = 'Soluble', degradability='Undegradable', organic=False)
    Nd2O3 = qs.Component.from_chemical('Nd2O3', chemical=Nd2O3, particle_size = 'Particulate', degradability='Undegradable', organic=False)
    Gypsum = qs.Component.from_chemical('Gypsum', chemical=Calcium_sulfate_dihydrate, particle_size = 'Particulate', degradability='Undegradable', organic=False)
    H2SO4 = qs.Component('H2SO4', search_ID='7664-93-9', particle_size = 'Soluble', degradability='Undegradable', organic=False)
    HNO3 = qs.Component('HNO3', search_ID='7697-37-2', particle_size = 'Soluble', degradability='Undegradable', organic=False)
    OA = qs.Component('OA', search_ID='144-62-7', particle_size = 'Soluble', degradability='Undegradable', organic=False) # anhydrous oxalic acid
    CO2 = qs.Component('CO2', search_ID='124-38-9', particle_size = 'Soluble', degradability='Undegradable', organic=False)
    CO = qs.Component('CO', search_ID='630-08-0', particle_size = 'Soluble', degradability='Undegradable', organic=False)
    NaOH = qs.Component('NaOH', search_ID='1310-73-2', particle_size = 'Soluble', degradability='Undegradable', organic=False)
    Na3PO4 = qs.Component.from_chemical('Na3PO4', chemical=Na3PO4, particle_size = 'Soluble', degradability='Undegradable', organic=False)
    uranium = qs.Component('U', search_ID='7440-61-1', particle_size = 'Soluble', degradability='Undegradable', organic=False)

    cmps = qs.Components([H2O, Hion, Neodymium, Nd2O3, Gypsum, H2SO4, HNO3, OA, CO2, CO, NaOH, Na3PO4, uranium]) # *cmps_default
    qs.set_thermo(cmps)

    return cmps
# cmps = create_components()
# print(cmps)

# cmps, Nd2O3 = create_components()
# print(Nd2O3.show(chemical_info=True))