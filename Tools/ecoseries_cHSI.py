import os, sys
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import arcpy
from arcpy import env
from arcpy.sa import *
try:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + "\\.site_packages\\riverpy\\")
    import config
    import fGlobal as fGl
except:
    print("ExceptionERROR: Missing RiverArchitect packages (riverpy).")

# Eco series analysis - cHSI
# This python script (1)

#########################
# User defined variables

case_name = "sfe_316"
scale_to_one = 0

probe_path_in = os.path.abspath("../SHArC/SHArea/" + case_name + "/probe.shp")
probe_path_out = os.path.abspath("../SHArC/SHArea/" + case_name + "/probe_cHSI.shp")
timeseries_path = "../00_Flows/" + case_name + "/flow_series_" + case_name + "_sample.xlsx"
figure_path = "../SHArC/SHArea/" + case_name + "/"

interptype = 'linear'
point_num = 3
colors = ["tab:blue", "tab:orange", "tab:green"]
#fish_periods = ["chju", "raju", "raad"]
fish_periods = ["chju"]
ind = 0

for fish_period in fish_periods:

    fish_name = fish_period[0:2]
    period = fish_period[2:4]

    if fish_name == 'ch':
        fish_full = 'Chinook Salmon'
    elif fish_name == 'ra':
        fish_full = 'Rainbow / Steelhead Trout'
    if period == 'ju':
        period_full = 'juvenile'
    elif period == 'ad':
        period_full = 'adult'

    fish_period_full = fish_full + ' - ' + period_full

    ######################
    # Reading probe points
    # Set environment settings
    env.workspace = os.path.abspath("../SHArC/SHArea/Rasters_" + case_name + "/no_cover")

    # Set local variables
    inPointFeatures = probe_path_in
    outPointFeatures = probe_path_out

    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    # Execute ExtractValuesToPoints
    rasters = arcpy.ListRasters("*", "tif")
    Flow = []
    cHSI = []
    ii = 0
    jj = 0

    for raster in rasters:
        if (raster[4:6] == fish_name) & (raster[6:8] == period):
            print(raster)
            Flow.append(fGl.read_Q_str(raster,'csi_'+fish_name+period))

            # get the coordinate system by describing a feature class
            dsc = arcpy.Describe(raster)
            coord_sys = dsc.spatialReference
            arcpy.management.DefineProjection(inPointFeatures, coord_sys)

            ExtractValuesToPoints(inPointFeatures, raster, outPointFeatures)

            rows = arcpy.UpdateCursor(outPointFeatures)

            for row in rows:
                # Fields from the table can be dynamically accessed from the row object.
                #   Here, the field is named targetField as I don't know your field name
                targetRow = row.RASTERVALU  # Assigns value of targetField to string
                #print(targetRow)
                if targetRow < 0:
                    targetRow = 0
                cHSI.append(targetRow)

            # Delete cursor and row objects to remove locks on the data
            #
            del row
            del rows

            ii += 1
            jj += 1

            arcpy.Delete_management(outPointFeatures)

    n = int(Flow.__len__())
    m = int(cHSI.__len__()/Flow.__len__())
    cHSI_orig = cHSI
    cHSI = cHSI_orig[slice(point_num-1, m*n, m)]

    Flow = np.append(Flow, [0])
    cHSI = np.append(cHSI, [0])

    Flow_new = np.linspace(np.min(Flow), np.max(Flow), num=10001, endpoint=True)
    f = interp1d(Flow, cHSI, kind=interptype)

    plt.figure(1)
    plt.plot(Flow, cHSI, marker="o", color=colors[ind], linewidth=0)
    plt.plot(Flow_new, f(Flow_new), color=colors[ind], label=fish_period_full)
    #plt.title(case_name + ', ' + fish_name + ', ' + period + ' at Point ' + str(point_num))
    plt.title(case_name+' at Point ' + str(point_num))
    plt.xlabel('Discharge ($m^3$/s)')
    plt.ylabel('cHSI at Point ' + str(point_num))
    if scale_to_one &(ind == fish_periods.__len__()-1):
        bottom, top = plt.ylim()
        plt.ylim(0, 1.3)
    plt.legend()
    plt.show()
    #plt.savefig(figure_path+case_name+'_'+ fish_name + period +'_cHSI_P' + str(point_num) + '_Q.svg')
    plt.savefig(figure_path + case_name + '_cHSI_P' + str(point_num) + '_Q.svg')
    plt.savefig(figure_path + case_name + '_cHSI_P' + str(point_num) + '_Q.pdf')

    #########################
    # Reading flow timeseries

    f3 = pd.read_excel(timeseries_path, index_col=None, usecols="A")[3:].values.tolist()
    f4 = pd.read_excel(timeseries_path, indox_col=None, usecols="B")[3:].values.tolist()

    Date = np.array(f3).transpose()[0]
    Flow_series = np.array(f4).transpose()[0]
    Eco_series = f(Flow_series)

    plt.figure(2)
    plt.plot(Flow_series, 'k')
    #plt.title(case_name + ', ' + fish_name + ', ' + period)
    plt.title(case_name)
    plt.xlabel('Time (days)')
    plt.ylabel('Discharge ($m^3$/s)')
    bottom, top = plt.ylim()
    plt.ylim(0, top)
    plt.show()
    #plt.savefig(figure_path+case_name+'_'+ fish_name + period +'_Q_time.svg')
    plt.savefig(figure_path + case_name + '_Q_time.svg')
    plt.savefig(figure_path + case_name + '_Q_time.pdf')

    plt.figure(3)
    plt.plot(Eco_series, label=fish_period_full)
    #plt.title(case_name + ', ' + fish_name + ', ' + period + ' at Point ' + str(point_num))
    plt.title(case_name + ' at Point ' + str(point_num))
    plt.xlabel('Time (days)')
    plt.ylabel('cHSI at Point ' + str(point_num))
    if (ind == fish_periods.__len__()-1):
        bottom, top = plt.ylim()
        plt.ylim(0, 1.3*top)
    plt.legend()
    plt.show()
    #plt.savefig(figure_path+case_name+'_'+ fish_name + period +'_cHSI_P' + str(point_num) + '_time.svg')
    plt.savefig(figure_path + case_name + '_cHSI_P' + str(point_num) + '_time.svg')
    plt.savefig(figure_path + case_name + '_cHSI_P' + str(point_num) + '_time.pdf')

    #########################
    # Sequence-average plot

    length = Eco_series.__len__()
    windows = range(365, length-1, 365)

    seq_min_series = []
    seq_avg_series = []
    seq_min_10 = []
    seq_min_90 = []
    seq_avg_10 = []
    seq_avg_90 = []

    for window in windows:
        for ii in range(0, length - window + 1):
            seq_min_series.append(np.min(Eco_series[ii:ii+window]))
            seq_avg_series.append(np.average(Eco_series[ii:ii+window]))
    #    seq_min_10.append(np.percentile(seq_min_series, 10))
    #    seq_min_90.append(np.percentile(seq_min_series, 90))
        seq_avg_10.append(np.percentile(seq_avg_series, 10))
        seq_avg_90.append(np.percentile(seq_avg_series, 90))


    plt.figure(4)
    #plt.plot(seq_min_10)
    #plt.plot(seq_min_90)
    plt.plot(seq_avg_10, colors[ind]) #, label='10 Percentile')
    plt.plot(seq_avg_90, colors[ind], label= fish_period_full)

    # Patches
    m = []
    for i in range(seq_avg_10.__len__()):
        m.append(i)
    x = np.hstack(([m], [m[::-1]]))
    y = np.hstack(([seq_avg_10], [seq_avg_90[::-1]]))
    patches = []
    polygon = Polygon(np.transpose(np.vstack((x,y))), True)
    patches.append(polygon)

    p = PatchCollection(patches, alpha=0.4, facecolors=colors[ind])
    #p.set_array(np.array(colors))
    plt.gca().add_collection(p)

    #plt.title(case_name + ', ' + fish_name + ', ' + period + ' at Point ' + str(point_num))
    plt.title(case_name + ' at Point ' + str(point_num))
    plt.xlabel('Length of sequence (year)')
    plt.ylabel('Sequence-averaged cHSI at Point ' + str(point_num))
    plt.xlim(0, 5)
    if scale_to_one &(ind == fish_periods.__len__()-1):
        bottom, top = plt.ylim()
        plt.ylim(0, 1.3)
    plt.legend()
    plt.show()
    #plt.savefig(figure_path + case_name + '_' + fish_name + period + '_cHSI_P' + str(point_num) + '_seq_avg.svg')
    plt.savefig(figure_path + case_name + '_cHSI_P' + str(point_num) + '_seq_avg.svg')
    plt.savefig(figure_path + case_name + '_cHSI_P' + str(point_num) + '_seq_avg.pdf')

    ind += 1