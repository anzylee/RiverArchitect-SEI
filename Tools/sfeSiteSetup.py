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

'''


'''

class model():
    # Initializer / Instance attributes
    def __init__(self, case, reference_site=2090, case_site_database='input/case_site_database.csv', site_area_database='input/All SFE LOI Characteristics with MAF.xlsx', scenario_path='input/scenarios/', sharea_path="../SHArC/SHArea/", fish_names=['Chinook Salmon', 'Rainbow / Steelhead Trout'], fish_stages=['fry','juvenile','spawn']):
        '''
        Parameters
        ----------
        case : Integer
            Catchment case number to label and find data.
        reference_site : Integer, optional
            DESCRIPTION. The default is 2090.
        case_site_database : String, optional
            DESCRIPTION. The default is 'input/case_site_database.csv'.
        site_area_database : String, optional
            DESCRIPTION. The default is 'input/All SFE LOI Characteristics with MAF.xlsx'.
        scenario_path : String, optional
            DESCRIPTION. The default is 'input/scenarios/'.
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
        # Initialize case and reference IDs and areas
        #cs_db = pd.read_csv(case_site_database) ##FINISH##
        self.case = case
        self.case_name = "sfe_" + '{0:d}'.format(case)
        self.case_site = 4370 ##FINISH##
        self.ref_site = reference_site
        self.case_area = 6.529 ##FINISH##
        self.ref_area = 26.04775 ##FINISH##
        
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
        This function calculates the area in the provided raster that meets the conditional expression indicated. The conditional expression can be either a minimum value or a maximum value.
        
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
            Contains information on bankfull area, flow exceedance, and normalized values of each of the previous.
    
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
            
                case_sharea_path = self.sha_path + self.case_name + "/" + self.case_name + "_sharea_" + fish_period + ".xlsx"
            
                ######################
                # Bankfull wetted area defined as the area with water present at highest flow in Flow
                # It is calculated for the first dataset and then remains constant for all following species/life stage datasets of the same river segment
                if ind == 0:
                    BfQ_area, dfSHArea_temp = self.ecoAreaDF(case_sharea_path=case_sharea_path, calc_BfQarea=True) # defined in sfeTools
                    
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
    
    def calculateSuccess(self, scenarios, eco_threshold=0.1, life_stage_period=False, verbose=False, df='2-d hydrodynamic'):
        '''

        Parameters
        ----------
        scenarios : List of integers
            List of scenarios to rank.
        eco_threshold : Float, optional
            A number between 0 and 1 that indicates the threshold to which the ecoseries habitat area to bankfull area ratio should be compared. The default is 0.1.
        life_stage_period : Dictionary of dictionaries, optional
            A dictionary that contains a dictionary of tuples for each species. Each species dictionary contains a tuple of length 3 for each life stage. The tuple indicates: ((start month, start day), (end month, end day), Including). Including is a boolean variable that indicates whether the date range to be evaluated is inside the date range indicated or outside the date range. The default for life_stage_period is False, which indicates the default dictionary should be used.
        verbose : Boolean, optional
            Indicates if progress of the analysis should be printed in the dialog. The default is False.

        Returns
        -------
        ecoseries_success : Pandas DataFrame
            Contains information on the number of successful days during each annual date range period when compared to the eco_threshold for all scenario, life stage, and fish name combinations.

        '''
        # Check to see if a life stage period dictionary is provided and if not, use the default
        if ~isinstance(life_stage_period, dict):
            life_stage_period = {}
            life_stage_period['Chinook Salmon'] = {}
            life_stage_period['Chinook Salmon']['fry'] = ((1,1),(12,31),True)
            life_stage_period['Chinook Salmon']['juvenile'] = ((1,1),(12,31), True)
            life_stage_period['Chinook Salmon']['spawn'] = ((1,1),(12,31), True)
            life_stage_period['Rainbow / Steelhead Trout'] = {}
            life_stage_period['Rainbow / Steelhead Trout']['fry'] = ((1,1),(12,31), True)
            life_stage_period['Rainbow / Steelhead Trout']['juvenile'] = ((6,16),(11,30), True)
            life_stage_period['Rainbow / Steelhead Trout']['spawn'] = ((1,1),(12,31), True)
        
        # Initialize DataFrame that will contain number of successes
        ecoseries_success = pd.DataFrame()
            
        for s in scenarios:
            
            if verbose: print(s)
            
            dfSHArea, dfTimeSeries, self.streamflowUnits, dfEcoseries = self.processData(s)
            
            years = np.unique(dfTimeSeries.index.year.values)
            
            for f in life_stage_period:
                
                if verbose: print(f)
                
                for d in life_stage_period[f]:
                    
                    if verbose: print(d)
                    
                    # Create list of dates to include for success ranking for each fish species and life stage combination
                    valid_dates = []
                    for y in years:
                        valid_dates.extend(list(pd.date_range(start=str(life_stage_period[f][d][0][0])+'/'+str(life_stage_period[f][d][0][1])+'/'+str(y), end=str(life_stage_period[f][d][1][0])+'/'+str(life_stage_period[f][d][1][1])+'/'+str(y))))
                    
                    dfEcoseries_temp = dfEcoseries[f + ' - ' + d].copy()
                    # Crop Eco Series DataFrame to given dates. The last tuple designates whether the fish period is between the two dates or outside of the two dates.
                    if life_stage_period[f][d][2]:
                        dfEcoseries_temp = dfEcoseries_temp[dfEcoseries_temp.index.isin(valid_dates)]
                    else:
                        dfEcoseries_temp = dfEcoseries_temp[~dfEcoseries_temp.index.isin(valid_dates)]
                    
                    # Determine if each entry meets the Eco Series Threshold
                    dfEcoseries_temp = dfEcoseries_temp > eco_threshold
                    
                    # Count all entries above the threshold in each year
                    numSuccess = pd.DataFrame(dfEcoseries_temp.resample('A').sum()).rename(columns={f + ' - ' + d: "Successes"})
                    numSuccess['Case'] = self.case
                    numSuccess['Scenario'] = s
                    numSuccess['Fish name'] = f
                    numSuccess['Life stage'] = d
                    numSuccess['Fish name - Life stage'] = f + ' - ' + d
                    
                    ecoseries_success = pd.concat([ecoseries_success, numSuccess])
        
        return ecoseries_success
            
    
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
