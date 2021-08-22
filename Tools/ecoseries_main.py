# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 15:37:10 2021

@author: MM
"""
import sfeSiteSetup as sfe
import seaborn as sns

sns.set_style("whitegrid")

case = 316
scenario = 5

# Initialization
sfe_case = sfe.model(case)
dfSHArea, dfTimeSeries, streamflow_name, dfEcoseries = sfe_case.processData(scenario)

# Plotting
sfe_case.caseBankfullQtoALinePlot(dfSHArea)
sfe_case.streamflowScatterPlot(dfTimeSeries)
sfe_case.scenarioHabitatSeriesLinePlot(dfEcoseries)
sfe_case.plot_multiSeqAvg(dfEcoseries)

# Scenario Ranking
ecoseries_success = sfe_case.calculateSuccess([16,17,18], verbose=True)

x = sfe_case.polyArea('depthraster.tif', limit=5, minlim=False)