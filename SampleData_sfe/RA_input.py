#############################################################################################
# RA_input.py prepares input for River Architect from tuflow simulation result files
#   Written by Anzy Lee, Postdoctoral Scholar, Utah State University
#   Date: 10/21/2020
#############################################################################################

import os
import shutil as sl
import numpy as np
import arcpy
import pandas as pd
import sys
import csv
import subprocess
import matplotlib.pyplot as plt

#############################################################################################
# 1 Input variables: case_name, create_folders, copy_depths_velocities, copy_terrain,
#                    [create_dmean, dmean], Q_all
case_name = 'sfe_2248' # 'VanillaC4', 'InphaseC4', 'OutphaseC4'
create_folders = 1
copy_depths_velocities = 1
copy_terrain = 1
[create_dmean, dmean] = [1, 0.001] # dmean = grain size in meters, 0.001m for sand bed

#############################################################################################
arcpy.env.overwriteOutput = True
tuflow_path = 'F:/tuflow_runs'

flow_stage = pd.read_excel(tuflow_path+'/'+case_name+'_tuflow/flow_stage_depth.xlsx')
Q_all = flow_stage.Q_cms # From tuflow_run

params = pd.read_excel('F:/tuflow-wf_python3/parameters_sfe.xlsx')
for ii in range(0,params.site.__len__()+1):
    if params.site[ii] == case_name:
        ind = ii
        break
n = params.n[ind]
S = params.Slope[ind]
grid_name = case_name + '.asc'

#############################################################################################

if create_folders == 1:
    os.mkdir('./01_Conditions/'+case_name)
    os.mkdir('./01_Conditions/' + case_name+'/depth')
    os.mkdir('./01_Conditions/' + case_name+'/velocity')

for ii in range(0,Q_all.__len__()): #np.array([9,19]): # range(0,Q_all.__len__())
    zero1 = ''
    zero2 = ''

    s = str(ii+1)
    s_len = s.__len__()

    for ind in range(0,3-s_len):
        zero1 = zero1 + '0'
    Q_num = zero1 + s

    Q_format = f'{Q_all[ii]:09.3f}'
    Q_format = Q_format.replace('.', '_')
    s_len = Q_format.__len__()

    for ind in range(0,10-s_len):
        zero2 = zero2 + '0'
    Q_value = zero2 + Q_format

    if copy_depths_velocities == 1:

        # depth: copy
        in_tif = os.path.abspath(tuflow_path+"/"+case_name+"_tuflow/results/clipped/T1_"+s+"_d_clip.tif")
        out_raster = os.path.abspath("./01_Conditions/"+case_name+"/depth/h"+Q_value+".tif")
        arcpy.CopyRaster_management(in_tif, out_raster)

        # velocity: flt to tif
        in_tif = os.path.abspath(tuflow_path+"/"+case_name+"_tuflow/results/clipped/T1_"+s+"_V_clip.tif")
        out_raster = os.path.abspath("./01_Conditions/"+case_name+"/velocity/u"+Q_value+".tif")
        arcpy.CopyRaster_management(in_tif, out_raster)

if copy_terrain == 1:
    in_ascii = os.path.abspath(tuflow_path+"/"+case_name+"_tuflow/model/grid/"+grid_name)
    out_raster = os.path.abspath("./01_Conditions/"+case_name+"/dem.tif")
    rasterType = "FLOAT"
    arcpy.ASCIIToRaster_conversion(in_ascii, out_raster, rasterType)

if create_dmean == 1:
    in_raster = os.path.abspath("./01_Conditions/"+case_name+"/dem.tif")
    out_raster = os.path.abspath("./01_Conditions/"+case_name+"/dmean0.tif")
    arcpy.CopyRaster_management(in_raster, out_raster)

    #arcpy.gp.CheckOutExtension("Spatial")
    #arcpy.gp.workspace = os.path.abspath("./01_Conditions/"+case_name)
    from arcpy.sa import *
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = os.path.abspath("./01_Conditions/"+case_name)
    in_raster = "dmean0.tif"
    out_raster = "dmean.tif"
    out_minus = Minus(in_raster, in_raster)
    out_dmean = Plus(out_minus,dmean)
    out_dmean.save(out_raster)
    arcpy.Delete_management(in_raster)
