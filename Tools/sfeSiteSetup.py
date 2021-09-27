# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 15:37:47 2021

@author: MM
"""
import fiona
import rasterio
import rasterio.features
from shapely.ops import transform
from shapely.geometry import shape, mapping
from shapely.geometry.multipolygon import MultiPolygon
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import seaborn as sns
import tkinter as tk
from datetime import timedelta
import sfeGUI

'''
To do:
Input dialog for initialization
- 2-d hydrodynamic, rct probe, or rct discharge
- Reference site, fish names, period, case site

Should the ecothreshold also vary with fish name and life stage?
    2-d hydrodynamic case: allow for default or table input

Graphing for RCT Tables will be by each combination, not combined as in the 2-D hydrodynamic case. Have a limit of 3 columns maximum.

ADD IN PNG OPTION FOR ALL GRAPHS

For the annual number of successful days, currently the year is January - December, excluding any days that are not in the range of days/months defined in the dictionary. We can change this to be the California water year, which would be October-September of the next year, with the same days excluded as before. Same with consecutive, etc YES, CALIFORNIA WATER YEAR, BAKED IN.

Greater than (>)

For the consecutive case with two thresholds, I am using the last date that the flow drops below the higher threshold, is that correct? THIS IS CORRECT, SHOW DECISIONS IN GUI

What should we do about maxflow for RCT option? MAKE IT AN OPTION, FOR RCT DO NOT LIMIT BY MAX FLOW, FOR 2-D HYDRODYNAMIC DO LIMIT BY MAX FLOW

Only two of the reference sites have streamflow from WEAP, what should we do? EMAIL CHUCK ABOUT RERUNNING WEAP WITH THE 7 REFERENCE STREAMS

For the discharge option, I am normalizing by the ratio of reference site to case site. Is that correct? DISCHARGE OPTION DOES NOT REQUIRE ANY CONVERSIONS (EXCEPT MAYBE FOR VISUALIZATION). PROBE REQUIRES CONVERSION OF STREAMFLOW FOR CALCULATING PROBE DEPTH, BUT DOES NOT REQUIRE CONVERSION OF THE RCT DEPTH THRESHOLD.

There are 5 points in the probe file. Does that mean you want me to check the depth at each of those points? Does it have to meet the threshold at all the points or just one of them? JUST CHECK THE MOST UPSTREAM ONE, PROBE.SHP IS UPDATED, update code to just use the first point

2-D HYDRODYNAMIC SHOULD ONLY HAVE ONE THRESHOLD OPTION (actually, add in a second window in the custom option that allows for flow or depth)

QUESTIONS:

It would be great to get the names of the watersheds/sites, so we could include them in the drop-down lists
'''

class model():
    # Initializer / Instance attributes
    def __init__(self, site_area_database='input/All SFE LOI Characteristics with MAF.xlsx', scenario_path='input/scenarios/', sharea_path="../SHArC/SHArea/", fish_names=['Chinook Salmon', 'Rainbow / Steelhead Trout'], fish_stages=['fry','juvenile','spawn']):
        '''
        Parameters
        ----------
        case : Integer
            Site of interest case number to label and find data.
        reference_site : Integer, optional
            Catchment ID for comparable catchment with empirical data available. The default is 2090.
        case_site_database : String, optional
            Database that contains information connecting the site case number to the catchment ID of the catchment in which it is located. The default is 'input/case_site_database.csv'.
        site_area_database : String, optional
            Database that contains contributing area data for each catchment ID. The default is 'input/All SFE LOI Characteristics with MAF.xlsx'.
        scenario_path : String, optional
            Location of streamflow timeseries for all scenarios to be analyzed. The default is 'input/scenarios/'.
        sharea_path : String, optional
            DESCRIPTION. The default is "../SHArC/SHArea/".
        fish_names : List of strings, optional
            Names for each of the species to be analyzed. The default is ['Chinook Salmon', 'Rainbow / Steelhead Trout'].
        fish_stages : List of strings, optional
            Life stages for fish to be analyzed. The default is ['fry','juvenile','spawn'].

        Returns
        -------
        Model object.
        '''
        # Set the analysis type
        self.modelType = sfeGUI.chooseModel()
        
        # Set the streamflow reference watershed
        self.ref_site = sfeGUI.chooseRefSite()
        
        # Set the hydrodynamic case site watershed
        if (self.modelType == 'rct probe') or (self.modelType == '2-d hydrodynamic'):
            self.case, self.case_site = sfeGUI.chooseCaseSite()
            if (self.modelType == 'rct probe'):
                self.case_name = "rct-probe_sfe-" + '{0:d}'.format(self.case)
            else:
                self.case_name = "2d-hydro_sfe-" + '{0:d}'.format(self.case)
        
        # Set the RCT paradigm model watershed
        if (self.modelType == 'rct probe') or (self.modelType == 'rct discharge'):
            rct_site = sfeGUI.chooseRCTParadigm()
            self.rct_site = rct_site
            
            if self.modelType == 'rct discharge':
                self.case_site = rct_site
                self.case_name = 'rct-discharge_' + '{0:d}'.format(self.case_site)
            else:
                self.case_name = self.case_name + '_{0:d}'.format(rct_site)
        
        self.life_stage_period = sfeGUI.setThresholds(self)
        
        # Initialize case and reference IDs and areas
        self.case_area = 6.529 ##FINISH##
        self.ref_area = 26.04775 ##FINISH##
        self.cfs = True
        
        # Initialize pathnames and headers
        self.s_path = scenario_path
        self.t_headers = list(pd.read_csv('input/SFER_Instance_Results_Template_long.csv').columns)
        self.sha_path = sharea_path
        self.fig_path = self.sha_path + self.case_name + "/"
        
        # Initialize species names and life stages
        self.f_name = fish_names
        self.f_stage = fish_stages
    
    def polyArea(self, raster_name, raster_location, limit=0, minlim=True, save_shp=False):
        '''
        Calculates the area in the provided raster that meets the conditional expression indicated. The conditional expression can be either a minimum value or a maximum value. There is also an option to save the shapefile of the calculated area.
        
        Parameters
        ----------
        raster_name : String
            Original raster to evaluate area.
        raster_location : String, optional
            File path to the where the raster that is being evaluated is stored.
        limit : float
            The value above or below which the area will be summed.
        minlim : boolean, optional
            Indicates whether the limit is a minimum vlaue (True) or a maximum value (False). The default is True.
        save_shp : boolean, optional
            If True, the shapefile associated with the polygon will be saved with the same name and path as the original raster. The default is False.
    
        Returns
        -------
        calculatedArea : float
            The area of the original raster that meets the limit requirement.
        '''
        filepath = os.path.abspath(raster_location)
        
        # Retrieve cell size and coordinate transformations for calculating area and generating shp
        tfw_metadata = np.loadtxt(os.path.join(filepath, raster_name + ".tfw"))
        cell_x = tfw_metadata[0]
        cell_y = tfw_metadata[3]
        top_left_x = tfw_metadata[4]
        top_left_y = tfw_metadata[5]
                   
        # Open raster and apply upscaling to counteract smoothing of the polygonize function below
        with rasterio.open(os.path.join(filepath, raster_name + ".tif")) as dem_src:
            upscale_factor = 2
            
            dtm_pre_arr = dem_src.read(
                out_shape=(
                    dem_src.count,
                    int(dem_src.height * cell_x * upscale_factor),
                    int(dem_src.width * cell_x * upscale_factor)
                )
            )
        
        # Determine if limit value is a minimum or maximum
        if minlim:
            outRas = dtm_pre_arr > limit
        else:
            outRas = dtm_pre_arr < limit
        shapes = list(rasterio.features.shapes(outRas.astype(int),mask=outRas,))
        
        # Create a MultiPolygon geometry with shapely.
        # shape() function is used to translate between GeoJSON and shapely geometry type
        polygons = [shape(geom) for geom, value in shapes
                    if value == 1]
        multipolygon = MultiPolygon(polygons)
        
        # Transform to crs from tfw file for accurate area and geolocation
        multipolygon2 = transform(lambda x, y, z=None: ((x/upscale_factor+top_left_x-cell_x/2), (-1*(y/upscale_factor+top_left_y+cell_y/2)+2*top_left_y)), multipolygon)
        
        # Write the record to an output shapefile with fiona if save_shp
        if save_shp:
            outPolygons = os.path.join(filepath, raster_name + ".shp")
            
            shp_schema = {
                'geometry': 'MultiPolygon',
                'properties': {'pixelvalue': 'int'}
            }
            
            with fiona.open(outPolygons, 'w', 'ESRI Shapefile', schema=shp_schema, crs=str(dem_src.crs)) as shp:
            # mapping() function is used to translate between shapely geometry type and GeoJSON dict 
                shp.write({
                    'geometry': mapping(multipolygon2),
                    'properties': {'pixelvalue': int(1)}
                })
        
        # Calculated area
        calculatedArea = np.float64(multipolygon2.area)
    
        return calculatedArea
    
    def ecoAreaDF(self, case_sharea_path, calc_BfQarea=True, BfQ_area=np.nan):
        '''
        Calculates the bankfull area, discharge normalized to bankfull discharge, and habitat area normalized to bankfull area given the case site and filename for an excel containing habitat area to discharge information. The excel should contain discharge (labeled as "Discharge") and calculated habitat areas (labeled as "Calculated Area") up to the maximum streamflow for the case site before water begins to leave the channel.
        
        Parameters
        ----------
        case_sharea_path : String
            File path to read in SHARrC data for this case.
        calc_BfQarea : Boolean
            Conditional to determine if BfQ_area should be calculated within the function or is given
        BfQ_area : float
            If calc_BfQarea is False, this is the value used to normalize the Calculated Area. The default is NaN.
        BfQ_depth_limit : float, optional
            The depth above which the area will be summed. The default is 0.
    
        Returns
        -------
        BfQ_area : float
            The area of the original depth raster that meets the limit requirement.
        dfSHArea_temp : pandas DataFrame
            Contains information on bankfull area, flow exceedance, and habitat area, and normalized values of each of the previous.
    
        '''
        # Read in SHARrC data
        dfSHArea_temp = pd.read_excel(case_sharea_path, index_col=None, header=0, skiprows=[0,2], usecols=[1,2,3,4,5])
        # Add an entry for no flow
        dfSHArea_temp = dfSHArea_temp.append(pd.DataFrame(data=[[0,'','',100,0]], index=[dfSHArea_temp.shape[0]], columns=list(dfSHArea_temp.columns)))
        dfSHArea_temp = dfSHArea_temp.rename(columns={'Calculated Area': 'Habitat area'})
        
        # Calculate bankfull area if not provided
        if calc_BfQarea:
            BfQ_depth = 'h' + f'{dfSHArea_temp.Discharge[0]:010.3f}'.replace('.', '_')
            BfQ_area = self.polyArea(raster_name=BfQ_depth, raster_location="../01_Conditions/"+ self.case_name, limit=0)
        
        # Non-dimensionalization
        dfSHArea_temp['Ratio of discharge to bankfull discharge'] = dfSHArea_temp['Discharge'].values / dfSHArea_temp['Discharge'][0]
        # Normalize all over Bank Full area
        dfSHArea_temp['Ratio of habitat area to bankfull area'] = dfSHArea_temp['Habitat area'].values / BfQ_area
        
        return BfQ_area, dfSHArea_temp
    
    def loadTimeSeries(self, timeseries_path, headers, siteid, maxflow, convert_to_cms=True, streamflow_normval=1):
        '''
        Loads a streamflow timeseries, extracts timeseries values that correspond only to the reference ID, normalizes by dividing by the provided factor, and converts from cubic feet per second to cubic meters per second if indicated. The source data should have a column that indicates the reference ID and the streamflow in cubic feet per second, without any headers. The headers must be provided in the function inputs and must label the streamflow as "Streamflow (cfs)" and the reference ID as "LOI".
        
        Parameters
        ----------
        timeseriespath : String
            Filepath for the streamflow time series.
        headers : List of strings
            Header names for the streamflow time series.
        siteid : Integer
            Site ID number.
        maxflow : Float
            Flow above which values should be set to 0.
        convert_to_cms : Boolean, optional
            Indication if the data should be converted from cubic feet per second to cubic meter per second. The default is True.
        streamflow_normval : Float, optional
            Value by which the streamflow timeseries should be normalized. The default is 1.
    
        Returns
        -------
        dfTimeSeries : Pandas DataFrame
            Streamflow timeseries with date information.
        '''
        # Import time series for current scenario
        dfTimeSeries = pd.read_csv(timeseries_path, header=None, names=self.t_headers, index_col='Date')
        dfTimeSeries.index = pd.to_datetime(dfTimeSeries.index, format=" %m/%d/%Y")
        
        # Select the subset of TimeSeries that corresponds to current Site_ID
        dfTimeSeries = dfTimeSeries[dfTimeSeries['LOI'] == siteid]
        dfTimeSeries = dfTimeSeries['Streamflow (cfs)']
        
        # Complete any conversion or normalization
        if convert_to_cms:
            self.cfs = False
            conversion = 0.0566
            units = '($m^3$/s)'
        else:
            conversion = 1
            units = '($f^3$/s)'
    
        dfTimeSeries = (np.floor(dfTimeSeries*1000)/1000*conversion) / streamflow_normval
            
        # Set any discharge above the maxflow (generally maxflow is set to bankfull discharge) to 0
        dfTimeSeries[(dfTimeSeries > maxflow)] = 0
        
        # Convert to DataFrame
        dfTimeSeries = dfTimeSeries.to_frame(name='Discharge '+units)
        
        return dfTimeSeries
    
    def ecoTimeSeries(self, dfTimeSeries, interp_func, series_name, eco_normval):
        '''
        
        
        Parameters
        ----------
        dfTimeSeries : Pandas DataFrame
            Streamflow timeseries with date information.
        interp_func : Function
            The interpolation function used to relate streamflow to habitat area, generally linear interpolation using scipy.
        series_name : String
            Name of the time series, generally in the format of 'Fish name - life stage'.
        eco_normval : Float
            Value by which the habitat area should be normalized, generally the bankfull area.
    
        Returns
        -------
        eco_series : Pandas DataFrame
            Normalized habitat area as a time series.
    
        '''
        eco_series = pd.DataFrame()
        temp_series = dfTimeSeries.apply(interp_func) / eco_normval
        eco_series[series_name] = np.squeeze(temp_series.values)
        eco_series.index = temp_series.index
        
        return eco_series
    
    def processData(self, scenario, verbose=False):
        '''
        This function processes the data for a given set of fish names and life stages for the scenario input.
        
        Parameters
        ----------
        scenario : Integer
            ID for the scenario streamflow time series evaluated.

        Returns
        -------
        dfSHArea : pandas DataFrame
            Contains information on bankfull area, flow exceedance, and normalized values of each of the previous.
        dfTimeSeries : Pandas DataFrame
            Streamflow timeseries with date information for the scenario selected.
        streamflow_name : String
            Either 'Streamflow (cfs)' or 'Streamflow (cms)', depending on whether the streamflow has been converted to cubic meters per second from cubic feet per second. Default is 'Streamflow (cms)'.
        dfEcoseries : Pandas DataFrame
            Normalized habitat area for the scenario selected as a time series.

        '''
        
        timeseries_path = self.s_path + "Scenario"+'{0:03d}'.format(scenario)+".csv"
        
        ind = 0
        dfSHArea = pd.DataFrame()
        dfEcoseries = pd.DataFrame()
        dfEcoseries_long = pd.DataFrame()
        
        for fname in self.f_name:
            for fstage in self.f_stage:
                fish_period = (fname[0:2] + fstage[0:2]).lower()
                
                if verbose: print(fish_period)
                
                # Excel spreadsheet that contains data relating streamflow (column labeled 'Discharge') to habitat area (column labeled 'Calculated Area') for each fish name and life stage combination. This excel spreadsheet is created by the user and can be created using RiverArchitect for proper formatting.
                case_sharea_path = self.sha_path + self.case_name + "/" + self.case_name + "_sharea_" + fish_period + ".xlsx"
            
                ######################
                # Bankfull wetted area defined as the area with water present at highest flow in Flow
                # It is calculated for the first dataset and then remains constant for all following species/life stage datasets of the same river segment
                if ind == 0:
                    BfQ_area, dfSHArea_temp = self.ecoAreaDF(case_sharea_path=case_sharea_path, calc_BfQarea=True)
                    
                    # Read in flow timeseries
                    dfTimeSeries = self.loadTimeSeries(timeseries_path=timeseries_path, headers=self.t_headers, siteid=self.ref_site, maxflow=dfSHArea_temp['Discharge'].values[0], streamflow_normval=self.ref_area/self.case_area)
                    self.streamflowUnits = dfTimeSeries.columns[0]
                    
                else:
                    BfQ_area, dfSHArea_temp = self.ecoAreaDF(case_sharea_path=case_sharea_path, calc_BfQarea=False, BfQ_area=BfQ_area)
                
                if verbose: print(BfQ_area)
                
                # Create a regression equation to relate Discharge to Habitat Area
                interpf = interp1d(dfSHArea_temp['Discharge'].values, dfSHArea_temp['Habitat area'].values, kind='linear')
                
                dfEcoseries_temp = self.ecoTimeSeries(dfTimeSeries, interp_func=interpf, series_name=fname + ' - ' + fstage, eco_normval=BfQ_area)
                dfEcoseries = pd.concat([dfEcoseries, dfEcoseries_temp], axis=1)
        
                # Add in specific information for each loop
                dfSHArea_temp['Case'] = self.case
                dfSHArea_temp['Scenario'] = scenario
                dfSHArea_temp['Fish name'] = fname
                dfSHArea_temp['Life stage'] = fstage
                dfSHArea_temp['Fish name - Life stage']= fname + ' - ' + fstage
                
                dfEcoseries_temp = dfEcoseries_temp.rename(columns={fname + ' - ' + fstage: 'Habitat area / Bankfull area'})
                dfEcoseries_temp['Case'] = self.case
                dfEcoseries_temp['Scenario'] = scenario
                dfEcoseries_temp['Fish name'] = fname
                dfEcoseries_temp['Life stage'] = fstage
                dfEcoseries_temp['Fish name - Life stage']= fname + ' - ' + fstage
                
                dfSHArea = pd.concat([dfSHArea,dfSHArea_temp])
                dfEcoseries_long = pd.concat([dfEcoseries_long,dfEcoseries_temp])
                ind += 1
        
        return dfSHArea, dfTimeSeries, self.streamflowUnits, dfEcoseries, dfEcoseries_long
        
    def calculateSuccess(self, scenarios, eco_threshold=0.1, verbose=False, analysis='2-d hydrodynamic'):
        '''

        Parameters
        ----------
        scenarios : List of integers
            List of scenarios to rank.
            
        life_stage_period : Dictionary of dictionaries
            A dictionary that contains a dictionary of tuples for each species. Each species dictionary contains a tuple of length 3 for each life stage. The tuple indicates: ((start month, start day), (end month, end day), T1, T2, Consecutive, Lower).
            T1 : Float
                A variable that indicates the primary Riffle Crest Thalweg (RCT) determined flow threshold for discharge for the reference watershed above which success is indicated.
            T2 : Float
                A variable that indicates the secondary Riffle Crest Thalweg (RCT) determined flow threshold for discharge for the reference watershed above which success is indicated. The default for T2 is False, which indicates it is not used.
            Consecutive : Boolean
                A variable that indicates how the success is evaluated according to the two T variables.
                'True' - If T1 and T2 are provided, indicates that the start of the period is the last date when the discharge drops below T1 and the end of the period is the first date when the discharge drops below T2. If only T1 is provided, indicates that the start of the period is the date provided and the end of the period is the first date when the discharge drops below T1. The number of successful days is the number of consecutive days counted.
                'False' - If T1 and T2 are provided, indicates that any days that are between T1 and T2 during the period between the start and end dates are counted as successful days. If only T1 is provided, indicates that any days that are above T1 during the period between the start and end dates are counted as successful days.
            Lower : Boolean
                A variable that indicates if the T2 variable is used or not.
                'True' - T1 and T2 are used.
                'False' - Only T1 is used.
                
        eco_threshold : Float, optional
            A number between 0 and 1 that indicates the threshold to which the ecoseries habitat area to bankfull area ratio should be compared. The default is 0.1.
            
        verbose : Boolean, optional
            Indicates if progress of the analysis should be printed in the dialog. The default is False.
            
        analysis : String, optional
            Indicates the type of analysis the user would like to complete. The default is '2-d hydrodynamic'

        Returns
        -------
        ecoseries_success : Pandas DataFrame
            Contains information on the number of successful days during each annual date range period when compared to the eco_threshold for all scenario, life stage, and fish name combinations.

        '''

        # Initialize DataFrame that will contain number of successes
        ecoseries_success = pd.DataFrame()
            
        for s in scenarios:
            
            if verbose: print(s)
            
            if analysis == '2-d hydrodynamic':
                dfSHArea, dfTimeSeries, self.streamflowUnits, dfEcoseries, dfEcoseries_long = self.processData(s)
            
            else:
                
                if (analysis == 'rct probe') or (analysis == 'rct discharge'):
                        
                    timeseries_path = self.s_path + "Scenario"+'{0:03d}'.format(s)+".csv"
                    
                    # Excel spreadsheet that contains data relating streamflow (column labeled 'Discharge') to habitat area (column labeled 'Calculated Area') for each fish name and life stage combination. This excel spreadsheet is created by the user and can be created using RiverArchitect for proper formatting. The path used for the RCT Table discharge needs is the first available, which is for 'chfr' 
                    case_sharea_path = self.sha_path + self.case_name + "/" + self.case_name + "_sharea_chfr.xlsx"
                    
                    # Bankfull wetted area (BfQ_area) is unused, therefore fill with dummy of 1
                    BfQ_area, dfSHArea_temp = self.ecoAreaDF(case_sharea_path=case_sharea_path, calc_BfQarea=False, BfQ_area=1)
                        
                    # Read in flow timeseries for reference site and convert to equivalent streamflow in case site
                    dfTimeSeries = self.loadTimeSeries(timeseries_path=timeseries_path, headers=self.t_headers, siteid=self.ref_site, maxflow=dfSHArea_temp['Discharge'].values[0], streamflow_normval=self.ref_area/self.case_area)
                    self.streamflowUnits = dfTimeSeries.columns[0]
                    
                    dfEcoseries = pd.DataFrame(index=dfTimeSeries.index)
                    
                    # Create an ecoseries dataframe using discharge data for all entries, this will be the Ecoseries DataFrame for the 'rct discharge option'
                    for f in self.life_stage_period:
                        
                        for l in self.life_stage_period[f]:
                            # Convert reference site thresholds to case site thresholds
                            self.life_stage_period[f][l][2] = self.life_stage_period[f][l][2]/(self.ref_area/self.case_area)
                            self.life_stage_period[f][l][3] = self.life_stage_period[f][l][3]/(self.ref_area/self.case_area)
                  
                            dfEcoseries[f + ' - ' + l] = dfTimeSeries.values.ravel()
                        
                    # if analysis == 'rct probe':
                    #     return 'rct probe'
                        
                else:
                    print('The variable analysis much be input as one of the following options: \'2-d hydrodynamic\', \'rct probe\', or \'rct discharge\'. Your data was not evaluated.')
                    return
            
            years = np.unique(dfTimeSeries.index.year.values)
            
            for f in self.life_stage_period:
                
                if verbose: print(f)
                
                for l in self.life_stage_period[f]:
                    
                    if verbose: print(l)
                    
                    # Create list of dates to include for success ranking for each fish species and life stage combination
                    valid_dates = []
                    dfEcorisk_temp = dfEcoseries[f + ' - ' + l].copy()
                    
                    # Check to see if consecutive
                    if self.life_stage_period[f][l][4]:
                        
                        # If consecutive, then determine if lower threshold is provided
                        if self.life_stage_period[f][l][5]:
                            
                            # If two thresholds, find the last date that the flow is greater than the first threshold and the first date that the flow is less than or equal to the second threshold
                            dfEcorisk_temp_gt1 = dfEcorisk_temp[dfEcorisk_temp > self.life_stage_period[f][l][2]]
                            dfEcorisk_temp_lt2 = dfEcorisk_temp[dfEcorisk_temp <= self.life_stage_period[f][l][3]]
                            
                            for y in years:
                                # Create a subset for each model year
                                dfDates_temp = dfEcorisk_temp_gt1.copy()
                                dfDates_temp = dfDates_temp[dfDates_temp.index.isin(list(pd.date_range(start='10/01/'+str(y), end='09/30/'+str(y+1))))]
                                
                                # Check to make sure there are values greater than the first threshold
                                if dfDates_temp.shape[0] > 0:
                                    # startdate is the last entry of each year that is less than or equal to the first threshold
                                    startdate = dfDates_temp.index[-1]
                                    
                                    dfDates_temp = dfEcorisk_temp_lt2.copy()
                                    dfDates_temp = dfDates_temp[dfDates_temp.index.isin(list(pd.date_range(start=startdate, end='09/30/'+str(y+1))))]
                                    
                                    # Check to make sure there are values less than the second threshold
                                    if dfDates_temp.shape[0] > 0:
                                        # enddate is the first entry following startdate that is less than or equal to the second threshold
                                        enddate = dfDates_temp.index[0]
                                    else:
                                        enddate = pd.to_datetime('09/30/'+str(y+1))
                                    
                                    valid_dates.extend(list(pd.date_range(start=startdate+timedelta(days=1), end=enddate-timedelta(days=1))))
                                    
                        else:
                            # If only one threshold, find the first date that the data series is above the threshold and the first date it falls below the threshold
                            dfEcorisk_temp_lt = dfEcorisk_temp[dfEcorisk_temp <= self.life_stage_period[f][l][2]]
                            dfEcorisk_temp_gt = dfEcorisk_temp[dfEcorisk_temp > self.life_stage_period[f][l][2]]
                            
                            for y in years:
                                # Create a subset for each model year
                                dfDates_temp = dfEcorisk_temp_gt.copy()
                                dfDates_temp = dfDates_temp[dfDates_temp.index.isin(list(pd.date_range(start='10/01/'+str(y), end='09/30/'+str(y+1))))]
                                
                                # Check to make sure there are values greater than the threshold
                                if dfDates_temp.shape[0] > 0:
                                    # startdate is the first entry of each year that is greater than the threshold
                                    startdate = dfDates_temp.index[0]
                                    
                                    dfDates_temp = dfEcorisk_temp_lt.copy()
                                    dfDates_temp = dfDates_temp[dfDates_temp.index.isin(list(pd.date_range(start=startdate, end='09/30/'+str(y+1))))]
                                    
                                    # Check to make sure there are values less than the threshold
                                    if dfDates_temp.shape[0] > 0:
                                        # enddate is the first entry of each year that is less than the threshold
                                        enddate = dfDates_temp.index[0]
                                    else:
                                        enddate = pd.to_datetime('09/30/'+str(y+1))
                                        
                                    valid_dates.extend(list(pd.date_range(start=startdate, end=enddate-timedelta(days=1))))
                    
                    else:
                    
                        # If not consecutive, create a valid date subset for each model year based on given start and end dates
                        for y in years:
                            startyear = y
                            endyear = y
                    
                            # If start or end month are before October, put start or end date in the following calendar year
                            if self.life_stage_period[f][l][0][0] < 10:
                                startyear += 1
                            if self.life_stage_period[f][l][1][0] < 10:
                                endyear += 1
                                
                            # Determine if lower threshold is provided
                            if self.life_stage_period[f][l][5]:
                                # If two thresholds, count all days between start and end dates that are less than or equal to the first threshold and greater than the second threshold
                                dfEcorisk_temp_btw = dfEcorisk_temp[(dfEcorisk_temp <= self.life_stage_period[f][l][2]) & (dfEcorisk_temp > self.life_stage_period[f][l][3])]
                                dfDates_temp = dfEcorisk_temp_btw.copy()
                                
                            else:
                                # If only one threshold, count all days between start and end dates that are above the threshold
                                dfEcorisk_temp_gt = dfEcorisk_temp[dfEcorisk_temp > self.life_stage_period[f][l][2]]
                                dfDates_temp = dfEcorisk_temp_gt.copy()
                            
                            # Create dataframe with values meeting the threshold(s) between the dates indicated in the dictionary provided
                            dfDates_temp = dfDates_temp[dfDates_temp.index.isin(list(pd.date_range(start=str(self.life_stage_period[f][l][0][0])+'/'+str(self.life_stage_period[f][l][0][1])+'/'+str(startyear), end=str(self.life_stage_period[f][l][1][0])+'/'+str(self.life_stage_period[f][l][1][1])+'/'+str(endyear))))]
                                
                            valid_dates.extend(list(dfDates_temp.index))
                    
                    dfEcorisk_temp = dfEcorisk_temp[dfEcorisk_temp.index.isin(valid_dates)]
                    
                    # Count all entries above the threshold in each year
                    numSuccess = pd.DataFrame(dfEcorisk_temp.resample('A-OCT').count()).rename(columns={f + ' - ' + l: "Successes"})
                    numSuccess['Case'] = self.case
                    numSuccess['Scenario'] = str(s)
                    numSuccess['Fish name'] = f
                    numSuccess['Life stage'] = l
                    numSuccess['Fish name - Life stage'] = f + ' - ' + l
                    
                    ecoseries_success = pd.concat([ecoseries_success, numSuccess])
        
        return ecoseries_success, dfEcoseries
            
    
    def caseBankfullQtoALinePlot(self, dfSHArea):
        # Line plot of normalized habitat area over bankfull area versus normalized discharge to bankfull discharge
        g = sns.relplot(data=dfSHArea, x="Ratio of discharge to bankfull discharge", y="Ratio of habitat area to bankfull area", hue="Fish name", style="Fish name", markers=True, col='Life stage', kind='line')
        g.savefig(self.fig_path + self.case_name + '_SHArea_Q-obj.svg')
        g.savefig(self.fig_path + self.case_name + '_SHArea_Q-obj.pdf')
        plt.close()
        
    def streamflowScatterPlot(self, dfTimeSeries):
        # Scatterplot of streamflow timeseries
        g = sns.scatterplot(data=dfTimeSeries, x=dfTimeSeries.index, y=self.streamflowUnits, linewidth=0, s=10, color='k')
        g.figure.savefig(self.fig_path + self.case_name + '_Q_time-obj.svg')
        g.figure.savefig(self.fig_path + self.case_name + '_Q_time-obj.pdf')
        plt.close()
    
    def scenarioHabitatSeriesLinePlot(self, dfEcoseries):
        # Line plot of normalized habitat area calculated based on discharge regression output over time
        g = sns.relplot(data=dfEcoseries, x=dfEcoseries.index, y='Habitat area / Bankfull area', hue='Fish name', row='Life stage', kind="line", legend=False, aspect=2)
        g.axes[0][0].legend(self.f_name)
        plt.tight_layout()
        g.savefig(self.fig_path + self.case_name + '_SHArea_time-obj.svg')
        g.savefig(self.fig_path + self.case_name + '_SHArea_time-obj.pdf')
        plt.close()

    def plot_seqAvg(self, df, label, ax, CI=0.8, window=365):
        '''
        This function plots a time series using the Sequence Average method. This method plots the confidence interval for the average value over windows of size b*n, where b is the base window size and n is 1 to the length of the time series divided by b minus 1. The average value is calculated as the rolling average at each window size.
        
        Parameters
        ----------
        df : Pandas DataFrame
            A single column DataFrame that contains data that will be plotted using the Sequence Averaged method.
        label : String
            Name of the time series, generally in the format of 'Fish name - life stage'.
        ax : Figure axes
            Axes object on which to plot the sequence-averaged graph.
        CI : Float, optional
            The size of the confidence interval that will be plotted. The default is 0.8.
        window : Integer, optional
            The base window size that will be used to carry out the Sequence Averaged method, in days. The default is 365.
    
        Returns
        -------
        None.
    
        '''
        years = df.resample(str(window)+'d').mean().shape[0]
        seqAvgCI = np.zeros((years-1,2))
        half_CI = (1 - CI)/2
        for y in range(1,years):
            seqAvgCI[y-1,0] = df.rolling(window*y).mean().quantile(half_CI)
            seqAvgCI[y-1,1] = df.rolling(window*y).mean().quantile(1 - half_CI)
        
        ax.plot(np.arange(1,years), seqAvgCI[:,0])
        ax.plot(np.arange(1,years), seqAvgCI[:,1], list(ax.get_lines())[-1].get_color())
        ax.fill_between(x=np.arange(1,years), y1=seqAvgCI[:,0], y2=seqAvgCI[:,1], alpha=0.4, label=label)
    
    def plot_multiSeqAvg(self, dfEcoseries):
        '''
        This function carries out the Sequence Average plot for all fish name and life stage combinations in the dfEcoseries DataFrame.

        Parameters
        ----------
        dfEcoseries : Pandas DataFrame
            A multi-column DataFrame that contains data that will be plotted using the Sequence Averaged method.

        Returns
        -------
        None.

        '''
        g, axs = plt.subplots(len(self.f_stage), sharex=True, sharey=True)
        for s, stage in enumerate(self.f_stage):
            for f in self.f_name:
                # Sequence average
                self.plot_seqAvg(df=dfEcoseries[f + ' - ' + stage], label=f, ax=axs[s], CI=0.8, window=365)
            axs[s].set_title(stage)
            axs[s].set_xlim(1)
        
        plt.legend(loc='best')
        g.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.ylabel('Sequence-averaged habitat area / Bankfull area', labelpad=10)
        plt.xlabel('Sequence window size (years)')
        plt.grid(False)
        g.savefig(self.fig_path + self.case_name + '_SHArea_seq_avg-obj.svg')
        g.savefig(self.fig_path + self.case_name + '_SHArea_seq_avg-obj.pdf')
        plt.close()

    def plot_ecorisk(self, ecoseries_success=list(np.arange(5,156)), annual_d_threshold=40, plt_d_per_yr=True, plt_success_d=True, plt_success_yr=True, ):
        '''
        Parameters
        ----------
        ecoseries_success : Pandas DataFrame or List of ints, optional
            If Pandas DataFrame, contains information for each model period year on how many days of success there are for each species and life stage combination evaluated. If list of scenario IDs, the ecoseries_success DataFrame is created with the default inputs of the calculateSuccess function for the specified scenarios. The default is a list of all the scenarios run in the South Fork Eel River basin through WEAP. To plot results using non-default inputs, calculateSuccess should be run with the specified inputs and the resulting ecoseries_success DataFrame should be provided here for plotting.
        annual_d_threshold : Integer or Dictionary, optional
            If integer, indicates the number of days above which a successful year is determined. If dictionary, it must include an integer entry for each fish name and life stage combination that indicates the number of days above which a successful year is determined. The default is 40, meaning the ecorisk threshold must be met at least 40 days in a given year to indicate a successful year.
        plt_d_per_yr : Boolean, optional
            Determines if the line graph of number of days for each model year is plotted for all scenario, fish name, and life stage combinations. The default is True.
        plt_success_d : Boolean, optional
            Determines if the line graph of cumulative number of successful days over the entire model period is plotted for all scenario, fish name, and life stage combinations. The default is True.
        plt_success_yr : Boolean, optional
            Determines if the line graph of number of successful years is plotted for all scenario, fish name, and life stage combinations. The default is True.

        Returns
        -------
        ecoseries_success : Pandas DataFrame
            The dataset provided as input or created using the calculateSuccess default inputs depending on user input.
        cumulative_days : Pandas DataFrame
            A dataset with the number of cumulative successful days as determined using the eco_threshold provided in the calculateSuccess function over the entire model period for each scenario, fish name, and life stage combination.
        cumulative_years : Pandas DataFrame
            A dataset with the number of cumulative successful years as determined according to annual_d_threshold over the entire model period for each scenario, fish name, and life stage combination.

        '''
        
        # Check to see if ecoseries_success is list of scenarios, if it is, then run calculateSuccess
        if isinstance(ecoseries_success, list):
            ecoseries_success = self.calculateSuccess(ecoseries_success, verbose=True)
        
        dfColumns = ['Scenario','Case','Fish name','Life stage','Fish name - Life stage']
        
        cumulative_days = ecoseries_success.groupby(dfColumns).sum().reset_index()
        # Check to see if annual_d_threshold is a dictionary, if not, apply same scalar year threshold to all
        if not(isinstance(annual_d_threshold, dict)):
            cumulative_years = ecoseries_success[ecoseries_success['Successes'] >= annual_d_threshold]
            
        # If annual_d_threshold is a dictionary, apply the threshold by fish name and life stage
        else:
            cumulative_years = pd.DataFrame()
            # Loop through fish names
            for f in annual_d_threshold:
                # Loop through life stages
                for d in annual_d_threshold[f]:
                    cumulative_years_temp = ecoseries_success[ecoseries_success['Fish name - Life stage'] == f + ' - ' + d]
                    cumulative_years_temp = cumulative_years_temp[cumulative_years_temp['Successes'] >= annual_d_threshold[f][d]]
                    cumulative_years = pd.concat([cumulative_years, cumulative_years_temp])
                    
        cumulative_years = cumulative_years.groupby(dfColumns).count().reset_index()
        cumulative_years_temp = cumulative_days.copy()
        cumulative_years_temp['Successes'] = 0
        cumulative_years = pd.concat([cumulative_years, cumulative_years_temp]).drop_duplicates(dfColumns, keep='first')
        
        if plt_d_per_yr:
            g = sns.relplot(data=ecoseries_success, x='Scenario', y='Successes', hue=ecoseries_success.index.year, col='Life stage', row='Fish name', kind="line", aspect=2, palette=sns.color_palette("viridis_r", as_cmap=True), facet_kws={'sharey':'none'})
            
            g.savefig(self.fig_path + self.case_name + '_ecorisk_days-per-year.svg')
            g.savefig(self.fig_path + self.case_name + '_ecorisk_days-per-year.pdf')
            plt.close()
        
        if plt_success_d:
            g = sns.relplot(data=cumulative_days, x='Scenario', y='Successes', hue='Fish name', col='Life stage', row='Fish name', kind="line", aspect=2, legend=False, facet_kws={'sharey':'none'})
            
            g.savefig(self.fig_path + self.case_name + '_ecorisk_cumulative-days.svg')
            g.savefig(self.fig_path + self.case_name + '_ecorisk_cumulative-days.pdf')
            plt.close()
            
        if plt_success_yr:
            g = sns.relplot(data=cumulative_years, x='Scenario', y='Successes', hue='Fish name', col='Life stage', row='Fish name', kind="line", aspect=2, legend=False, facet_kws={'sharey':'none'})
            
            g.savefig(self.fig_path + self.case_name + '_ecorisk_success-years.svg')
            g.savefig(self.fig_path + self.case_name + '_ecorisk_success-years.pdf')
            plt.close()
        
        return ecoseries_success, cumulative_days, cumulative_years