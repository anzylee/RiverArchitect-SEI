import os, sys
import pandas as pd
import numpy as np
import simpledbf
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import arcpy
from arcpy import env
from arcpy.sa import *

try:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + "\\.site_packages\\riverpy\\")
    import config
    import fGlobal as fGl
except:
    print("ExceptionERROR: Missing RiverArchitect packages (riverpy).")

# Eco risk assessment - SHArea

#########################
# User defined variables

case_name = "sfe_316" # sfe_322, sfe_25
Site_ID = 2090
Cont_Area = 26.04775*2.59e6
#case_Area = 7901.08
Site_Cont_Area = 6.529*2.59e6
From_year = 1996
To_year = 2017
Num_success_criteria = 40

# Specify fish and eco-period
fish_periods = ["chju"] # ch - Chinook salmon, ju - juvenile
Eco_threshold = 0.1 #

baseline_path = "F:/sfe_hydrology/Eco_performance/FOR_WEAP/FOR_WEAP_SFER_110_v1.xlsx"
From_scenario = 16
To_scenario = 25
include_base = 0
if include_base == 1:
    scenario_num = ['base']
else:
    scenario_num = []
for ii in range(From_scenario, To_scenario+1):
    zero1 = ''
    s = str(ii + 1)
    s_len = s.__len__()

    for ind in range(0, 3 - s_len):
        zero1 = zero1 + '0'
    ii_str = zero1 + s
    scenario_num.append(ii_str)

#scenario_num = ['016', '026', '036', '046', '056', '066', '076', '086',
#                '096', '106', '116', '126', '136', '146']
#scenario_num = [ '134', '135', '136', '137', '138', '139', '140',
#                '141', '142', '143', '144', '145', '146', '147', '148', '149']
#scenario_num = ['016', '017', '018', '019', '020', '021', '022',
#                '023', '024', '025']
#scenario_num = [ '026', '027', '028', '029', '030', '031', '032',
#                '033', '034', '035']

SFE_LOI_path = "F:/sfe_hydrology/Eco_performance/All SFE LOI Characteristics with MAF.xlsx"

figure_path = "../SHArC/SHArea/" + case_name + "/"

interptype = 'linear'

#########################

ind = 0
colors = ["tab:blue", "tab:orange", "tab:green"]
Flow_Num = []
Cuml_Num_success = []
Num_success_year = []
Num_success_all = []

for fish_period in fish_periods:

    fish_name = fish_period[0:2]
    period = fish_period[2:4]

    if fish_name == 'ch':
        fish_full = 'Chinook Salmon'
        if period == 'ju':
            period_full = 'juvenile'
            #period_start = datetime.datetime(2010, 6, 16, 0, 0)
            #period_end = datetime.datetime(2010, 11, 29, 0, 0)
            period_start = datetime.datetime(2010, 1, 1, 0, 0)
            period_end = datetime.datetime(2010, 12, 31, 0, 0)
        elif period == 'ad':
            period_full = 'adult'
            period_start = datetime.datetime(2010, 1, 1, 0, 0)
            period_end = datetime.datetime(2010, 12, 31, 0, 0)
    elif fish_name == 'ra':
        fish_full = 'Rainbow / Steelhead Trout'
        if period == 'ju':
            period_full = 'juvenile'
            period_start = datetime.datetime(2010, 6, 16, 0, 0)
            period_end = datetime.datetime(2010, 11, 29, 0, 0)
        elif period == 'ad':
            period_full = 'adult'
            period_start = datetime.datetime(2010, 1, 1, 0, 0)
            period_end = datetime.datetime(2010, 12, 31, 0, 0)

    fish_period_full = fish_full + ' - ' + period_full

    sharea_path = "../SHArC/SHArea/" + case_name + "/" + case_name + "_sharea_" + fish_name + period + ".xlsx"

    ######################
    # Reading SHARrC data

    f1 = pd.read_excel(sharea_path, index_col=None, header=None,usecols="B")[3:].values.tolist()
    f2 = pd.read_excel(sharea_path, index_col=None, header=None,usecols="F")[3:].values.tolist()

    Flow = np.array(f1).transpose()[0]
    CalArea = np.array(f2).transpose()[0]

    Flow = np.append(Flow, [0])
    CalArea = np.append(CalArea, [0])

    ######################
    # Bankfull wetted area
    env.workspace = os.path.abspath("../SHArC/HSI/" + case_name)
    BfQ_hsi = "dsi_" + fish_period + fGl.write_Q_str(Flow[0]) + ".tif"

    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    # Execute ExtractValuesToPoints
    rasters = arcpy.ListRasters("*", "tif")

    for raster in rasters:
        if raster == BfQ_hsi:
            print(raster)

            outRas = Raster(BfQ_hsi) > -1

            outPolygons = "BfQ_polygon.shp"
            arcpy.RasterToPolygon_conversion(outRas, outPolygons)

            # Set local variables
            inZoneData = outPolygons
            zoneField = "id"
            inClassData = outPolygons
            classField = "id"
            outTable = "BfQ_polygon_table.dbf"
            processingCellSize = 0.01

            # Execute TabulateArea
            TabulateArea(inZoneData, zoneField, inClassData, classField, outTable,
                         processingCellSize, "CLASSES_AS_ROWS")

            BfQ_area_dbf = simpledbf.Dbf5(env.workspace + '\\' + outTable)
            BfQ_partial_area = BfQ_area_dbf.to_dataframe()
            BfQ_area = np.sum(np.array(BfQ_partial_area['Area']))

            del BfQ_area_dbf
            del BfQ_partial_area
            #del BfQ_area

            arcpy.Delete_management(outPolygons)
            arcpy.Delete_management(outTable)

    # Reverse
    #Flow = Flow[::-1]
    #CalArea = CalArea[::-1]

    # Non-dimensionalization
    print(BfQ_area)
    Norm_Flow = Flow / Flow[0]
    Norm_CalArea = CalArea / BfQ_area
    #os.system("pause")

    ######################

    Norm_Flow_new = np.linspace(np.min(Norm_Flow), np.max(Norm_Flow), num=10001, endpoint=True)
    Norm_f = interp1d(Norm_Flow, Norm_CalArea, kind=interptype)
    f = interp1d(Flow, CalArea, kind=interptype)

    plt.figure(1)
    plt.plot(Norm_Flow, Norm_CalArea, marker="o", color=colors[ind], linewidth=0)
    plt.plot(Norm_Flow_new, Norm_f(Norm_Flow_new), color=colors[ind], label= fish_period_full)
    #plt.title(case_name + ', ' + fish_name + ', ' + period)
    plt.title(case_name)
    plt.xlabel('Ratio of discharge to bankfull discharge')
    plt.ylabel('Habitat area / Bankfull area')

    if ind == fish_periods.__len__()-1:
        bottom, top = plt.ylim()
        if top > 0.5:
            plt.ylim(0, 1)
        else:
            plt.ylim(0, 0.5)
    plt.legend()
    plt.show()

    #########################
    # Reading flow timeseries

    # f3 = pd.read_excel(timeseries_path, index_col=None, usecols="A")[3:].values.tolist() #date time
    # f4 = pd.read_excel(timeseries_path, indox_col=None, usecols="B")[3:].values.tolist()
    #flows = ['base']
    flows = scenario_num
    ind_flow = flows.__len__()

    for flow in flows:
        Date = []

        if flow[0] == 'b':
            timeseries_path = baseline_path
            print('Reading excel file ...')
            xls = pd.ExcelFile(timeseries_path)
            f3 = pd.read_excel(xls, 'RO_CFS')
            f3.rename({str(Site_ID): 'LOI'}, axis=1, inplace=True)
            Flow_series = f3.LOI
            Flow_series = np.floor(Flow_series * 1000) / 1000 * 0.0566  # cfs to cms
            Flow_series = Flow_series / Cont_Area * Site_Cont_Area  # Area
            Flow_series[Flow_series > Flow[0]] = 0
            Date_str = f3.DTTM
            for tmp in Date_str:
                [a0,b0,c0] = [tmp.year, tmp.month, tmp.day]
                tmp1 = datetime.datetime(int(a0),  # year
                                                  int(b0),  # month
                                                  int(c0),  # day
                                                  0, 0)
                Date.append(tmp1)
        else:
            timeseries_path = "F:/sfe_hydrology/Eco_performance/SFER/09_30_2020_Results/Scenario"+flow+".csv"
            f3 = pd.read_csv(timeseries_path)
            Site_OX = np.array(f3).transpose()[3] == Site_ID
            Date_str = np.array(f3).transpose()[2][Site_OX]
            for tmp in Date_str:
                tmp = tmp[1:]
                #tmp1 = datetime.datetime.strptime(tmp, '%m/%d/%y')
                [a0, b0, c0] = tmp.split('/')
                tmp1 = datetime.datetime(int(c0),  # year
                                                  int(a0),  # month
                                                  int(b0),  # day
                                                  0, 0)
                Date.append(tmp1)

            Flow_series = np.array(np.array(f3).transpose()[5][Site_OX], dtype=np.float64)
            Flow_series = np.floor(Flow_series * 1000) / 1000 * 0.0566  # cfs to cms
            Flow_series = Flow_series / Cont_Area * Site_Cont_Area  # Area
            Flow_series[Flow_series > Flow[0]] = 0

        Year = np.arange(Date[0].year, Date[-1].year+1)
        Num_success = []
        Num_success_tmp = int(0)
        ind_year = 0
        ind_date = 0

        Eco_series = f(Flow_series)
        Norm_Eco_series = Eco_series / BfQ_area
        Norm_Flow_series = Flow_series / Flow[0]
        pre_year = Date[0].year

        for date in Date:
            year = date.year
            # print(date)
            period_start0 = datetime.datetime(year,  # year
                                              period_start.month,  # month
                                              period_start.day,  # day
                                              0, 0)
            period_end0 = datetime.datetime(year, # year
                                            period_end.month, # month
                                            period_end.day, # day
                                            0, 0)
            if period_end0 < period_start0:
                period_end0 = datetime.datetime(year+1,  # year
                                                period_end.month,  # month
                                                period_end.day,  # day
                                                0, 0)
            if (date >= period_start0) & (date <= period_end0):
                if Norm_Eco_series[ind_date] >= Eco_threshold:
                    Num_success_tmp += 1
            ind_date += 1

            if (pre_year < year) | (date == Date[-1]):
                Num_success.append(Num_success_tmp)
                Num_success_tmp = 0
                ind_year += 1

            pre_year = year

        Num_success = np.array(Num_success)
        Num_success_all = np.append(Num_success_all, Num_success)


        plt.figure(2)
        if flow[0] == 'b':
            plt.plot(Date, Flow_series, 'k-')
        else:
            plt.plot(Date, Flow_series)
        #plt.title(case_name + ', ' + fish_name + ', ' + period)
        plt.title(case_name)
        plt.xlabel('Time (year)')
        plt.ylabel('Discharge ($m^3$/s)')
        bottom, top = plt.ylim()
        plt.ylim(0, top)
        plt.show()
        plt.savefig(figure_path + case_name + '_ecorisk_flow_series.svg')
        plt.savefig(figure_path + case_name + '_ecorisk_flow_series.pdf')

        """
        plt.figure(3)
        if flow[0] == 'b':
            plt.plot(Date, Norm_Eco_series, 'k-', label='Scenario = '+flow)
        else:
            plt.plot(Date, Norm_Eco_series, label='Scenario = '+flow)
        #plt.title(case_name + ', ' + fish_name + ', ' + period)
        plt.title(case_name+', '+fish_period_full)
        plt.xlabel('Time (year)')
        plt.ylabel('Habitat area / Bankfull area')
        if ind == fish_periods.__len__()-1:
            bottom, top = plt.ylim()
            if top > 0.5:
                plt.ylim(0, 1)
            else:
                plt.ylim(0, 0.5)
        plt.legend()
        plt.show()
        
        plt.figure(4)
        #plt.bar(Year, Num_success, label= fish_period_full)
        if flow[0] == 'b':
            plt.plot(Year, Num_success, 'k-+',label='Scenario = '+flow)
        else:
            plt.plot(Year, Num_success, '-+', label='Scenario = ' + flow)
        plt.title(case_name+', '+fish_period_full)
        plt.xlabel('Year')
        plt.ylabel('Number of success days')
        if ind == fish_periods.__len__()-1:
            bottom, top = plt.ylim()
            plt.ylim(0, np.round(1.3*top))
        plt.legend()
        plt.show()
        plt.savefig(figure_path + case_name + '_num_days.svg')
        plt.savefig(figure_path + case_name + '_num_days.pdf')
        """
        Cuml_Num_success_tmp = 0
        Num_success_year_tmp = 0
        ind_num = 0
        for tmp_year in Year:
            if (tmp_year >= From_year) & (tmp_year <= To_year):
                Cuml_Num_success_tmp = Cuml_Num_success_tmp + Num_success[ind_num]
                if Num_success[ind_num] >= Num_success_criteria:
                    Num_success_year_tmp += 1
            ind_num += 1
        print('Scenario = '+flow)
        print(Cuml_Num_success_tmp)

        #Flow_Num.append([flow, Cuml_Num_success])
        Cuml_Num_success.append(Cuml_Num_success_tmp)
        Num_success_year.append(Num_success_year_tmp)

Cuml_Num_success = np.array(Cuml_Num_success)
argsort_Cuml = np.argsort(-(Cuml_Num_success)) # minus sign for reverse
Cuml_Num_success.sort()
sort_flows = np.array(flows)[argsort_Cuml]
Cuml_Num_success = Cuml_Num_success[::-1] # reverse

plt.figure(5)
plt.plot(sort_flows, Cuml_Num_success)
plt.xlabel('Ranked instream flow - water management scenarios')
plt.ylabel('Successful days from '+str(From_year)+' to '+str(To_year))
#plt.plot(flows, Cuml_Num_success_orig)


plt.figure(6)
Num_success_all = np.resize(Num_success_all, (ind_flow, ind_year))
Num_success_all = Num_success_all.transpose()
for ind in range(0, ind_year):
    plt.plot(sort_flows, Num_success_all[ind], label=str(Year[ind]))
plt.xlabel('Ranked instream flow - water management scenarios')
plt.ylabel('Successful days')
plt.legend(bbox_to_anchor=(1.05, 1.1), loc='upper left')
plt.tight_layout()
plt.show()

Num_success_year = np.array(Num_success_year)
argsort_Num = np.argsort(-(Num_success_year)) # minus sign for reverse
Num_success_year.sort()
sort_flows = np.array(flows)[argsort_Num]
Num_success_year = Num_success_year[::-1] # reverse

plt.figure(7)
plt.plot(sort_flows, Num_success_year)
plt.xlabel('Ranked instream flow - water management scenarios')
plt.ylabel('Successful years from '+str(From_year)+' to '+str(To_year))