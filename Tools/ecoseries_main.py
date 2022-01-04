# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 15:37:10 2021

@author: MM
"""
import sfeSiteSetup as sfe

# Initialization
sfe_case = sfe.model()

print(sfe_case.modelType)
print(sfe_case.ref_site)
print(sfe_case.case_name)
print(sfe_case.case_site)

# Ecorisk analysis and plotting to compare multiple scenarios
if sfe_case.compare_scenario:
    if sfe_case.modelType == '2-d hydrodynamic':
        ecoseries_success, cummulative_days, cummulative_years = sfe_case.plot_ecorisk_2dh(sfe_case.scenarios)

    elif sfe_case.modelType == 'rct probe':
        ecoseries_success, cumulative_days, cumulative_years = sfe_case.plot_ecorisk_rct(sfe_case.scenarios)

    elif sfe_case.modelType == 'rct discharge':
        ecoseries_success, cumulative_days, cumulative_years = sfe_case.plot_ecorisk_rct(sfe_case.scenarios)

# Plotting to visualize a single scenario
else:
    if sfe_case.modelType == '2-d hydrodynamic':
        # Analyze single scenario
        dfSHArea, dfTimeSeries, dfEcoseries, dfEcoseries_long = sfe_case.processData_2dh(sfe_case.scenarios)
        
        # Plotting
        sfe_case.caseBankfullQtoALinePlot(dfSHArea)
        sfe_case.streamflowScatterPlot(dfTimeSeries)
        sfe_case.habitatSeriesLinePlot(dfEcoseries_long)
        sfe_case.plot_multiSeqAvg_2dh(dfEcoseries, window=365, CI=0.8)
        
    elif sfe_case.modelType == 'rct probe':
        # Analyze single scenario
        dfHSI, dfTimeSeries, dfEcoseries, dfEcoseries_long = sfe_case.processData_rctprobe(sfe_case.scenarios)
        
        # Plotting
        sfe_case.QtoProbeDepthLinePlot(dfHSI)
        sfe_case.streamflowScatterPlot(dfTimeSeries)
        sfe_case.probeSeriesLinePlot(dfEcoseries)
        sfe_case.plot_multiSeqAvg_rct(dfEcoseries, window=365, CI=0.8)
        
    elif sfe_case.modelType == 'rct discharge':
        # Analyze single scenario
        dfTimeSeries, dfEcoseries = sfe_case.processData_rctdischarge(sfe_case.scenarios)
       
        # Plotting
        sfe_case.streamflowScatterPlot(dfTimeSeries)
        sfe_case.plot_multiSeqAvg_rct(dfEcoseries, window=365, CI=0.8)
    