# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 15:37:10 2021

@author: MM
"""
import sfeSiteSetup as sfe
import pandas as pd
import seaborn as sns
import numpy as np

# Initialization
scenario = 5
sfe_case = sfe.model()

print(sfe_case.modelType)
print(sfe_case.ref_site)
print(sfe_case.case_name)
print(sfe_case.case_site)

if sfe_case.modelType == '2-d hydrodynamic':
    dfSHArea, dfTimeSeries, dfEcoseries, dfEcoseries_long = sfe_case.processData_2dh(scenario)
    # Plotting
    sfe_case.caseBankfullQtoALinePlot(dfSHArea)
    sfe_case.streamflowScatterPlot(dfTimeSeries)
    sfe_case.scenarioHabitatSeriesLinePlot(dfEcoseries_long)
    sfe_case.plot_multiSeqAvg_2dh(dfEcoseries)
    ecoseries_success, cummulative_days, cummulative_years = sfe_case.plot_ecorisk_2dh(ecoseries_success=[5,6,7,8])
    
elif sfe_case.modelType == 'rct probe':
    dfTimeSeries, dfEcoseries, dfEcoseries_long = sfe_case.processData_rctprobe(scenario)
    # Plotting
    sfe_case.streamflowScatterPlot(dfTimeSeries)
    ecoseries_success, cumulative_days, cumulative_years = sfe_case.plot_ecorisk_rct(ecoseries_success=[5,6,7,8], plt_d_per_yr=True, plt_success_d=False, plt_success_yr=False)
    
else:
    dfTimeSeries, dfEcoseries = sfe_case.processData_rctdischarge(scenario)
    # Plotting
    sfe_case.streamflowScatterPlot(dfTimeSeries)
    