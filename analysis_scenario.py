import numpy as np 
import scipy as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import qsdsan as qs
import biosteam as bst
import os

from systems import create_system

# Scenario Analysis
def analysis_scenario(fununit, num_ind_REEs, figures_path):
    x = np.linspace(0.02/100, 0.9/100, 20) # wt fraction. REE content
    y = np.linspace(0.1*1e6, 2*1e6, 20) # M kg/hr. capacity

    resultCapacity = []
    resultREEContent = []
    resultMSP = []
    for i in x:
        for j in y:
            # dataCapacity, dataREEContent, dataMSP = run_analysis(fununit='PG', feedPG= j, REEcontent= i, num_ind_REEs=9,
            #     report='no', num_samples=1, uncertainty='no', sensitivity='no', parameter='technological', optimization='no', desire_target='no')
            
            
            sys, lca, tea= create_system(fununit=fununit, feedPG=j, REEcontent=i, num_ind_REEs=num_ind_REEs)
            flowsheet = qs.Flowsheet.flowsheet.default
            fs_stream = flowsheet.stream

            dataCapacity = j
            dataREEContent = i
            dataMSP = tea.solve_price(fs_stream.Ln2O3)
            
            resultCapacity.append(dataCapacity)
            resultREEContent.append(dataREEContent)
            resultMSP.append(dataMSP)

    X, Y = np.meshgrid(x*100, y/1e6)
    Z = np.array(resultMSP).reshape(len(y),len(y)).T

    plt.style.use('default')
    aspect_ratio_LtoW = 2880/3840
    cm_to_in = 1/2.54  # centimeters in inches
    width_one_col = 8.3 # cm. Width for a one column figure
    width_two_col = 17.1 # cm. Width for a two column figure
    max_length = 23.3 # cm. The maximum lenght a figure can be
    fig, ax = plt.subplots(figsize=(width_one_col*cm_to_in, width_one_col*aspect_ratio_LtoW*cm_to_in))
    cp = plt.contourf(X,Y,Z, levels=np.linspace(0,400,30), cmap='inferno_r', extend='max')
    fig.colorbar(cp, label=f'MSP ($/kg REO)', ticks=np.arange(0,370,80)) 

    cp2 = plt.contour(X,Y,Z, levels=[51.5, 51.5*2], colors='black', linestyles='dashed')
    # fmt = {}
    # strs = ['  Selling Price  ', '  2x Selling Price  '] # ' MSP  51.5 ', ' MSP  103 '
    # for l, s in zip(cp2.levels, strs):
    #     fmt[l] = s
    # ax.clabel(cp2, inline=True, fmt=fmt)

    ax.set(xlabel='REE Content of the PG (wt %)', ylabel='Capacity (M kg PG/hr)')
    fig.tight_layout()
    fig.savefig(os.path.join(figures_path, f'Scenario_analysis.tiff'), dpi=600)
