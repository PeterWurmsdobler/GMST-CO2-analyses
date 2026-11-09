# Global Holocene Mean Surface Temperature Reconstructions
#----------------------------------------------------------------------- 
#                World Data Service for Paleoclimatology, Boulder 
#                                  and 
#                     NOAA Paleoclimatology Program 
#             National Centers for Environmental Information (NCEI)
#----------------------------------------------------------------------- 
# Template Version 3.0
# Encoding: UTF-8
# NOTE: Please cite Publication, and Online_Resource and date accessed when using these data. 
# If there is no publication information, please cite Investigators, Title, and Online_Resource and date accessed. 
#
# Online_Resource: https://www.ncdc.noaa.gov/paleo/study/29712
#       Description: NOAA Landing Page
# Online_Resource: https://www1.ncdc.noaa.gov/pub/data/paleo/reconstructions/kaufman2020/readme-kaufman2020gmst.txt 
#       Description: NOAA location of the template
#
# Original_Source_URL:  
#       Description:
# 
# Description/Documentation lines begin with #
# Data lines have no #
#
# Data Type: Climate Reconstructions
#
# Dataset DOI: 10.25921/vzys-1280 
#
# Parameter_Keywords: air temperature
#--------------------
# Contribution_Date
#	Date: 2020-06-12
#--------------------
# File_Last_Modified_Date
#	Date: 2020-06-12
#--------------------
# Title
#	Study_Name: Global Holocene Mean Surface Temperature Reconstructions
#--------------------
# Investigators
#	Investigators: Kaufman, D.; McKay, N.; Routson, C.; Erb, M.; Dätwyler, C.; Sommer, P.; Heiri, O.; Davis, B. 
#--------------------
# Description_Notes_and_Keywords
#	Description: Global mean surface temperature (GMST) and associated uncertainties for the past 12,000 years, 
#	using five different statistical methods. 
#	The reconstructions are based on v1.0 of the Temperature 12k Database, available at: 
#	https://www.ncdc.noaa.gov/paleo/study/27330
#	
#	Data have been stored both as netCDF and csv files.  
#	These files are described below, and contain data for the following five methods: 
#	
#	SCC: Standard calibrated composite 
#	DCC: Dynamic calibrated composite 
#	GAM: Generalized additive model 
#	CPS: Composite plus scale 
#	PAI: Pairwise comparison (PaiCo) 
#	
#	FILES:
#	temp12k_alldata.nc: 
#	This is a netCDF file containing reconstruction data for all methods, both for global and latitude-zone means, 
#	stored in the following variables. 
#	age: Ages for the 100-year bins spanning from 0 to 12000 years before present (yr BP). "Present" is 1950 CE.  
#	[method]_globalmean: Global-mean temperature reconstructions for all 500 ensemble members of a given method.  
#	Format: Matrix of size 121 (ages) x 500 (ensemble members)  
#	[method]_latbands: Latitude-band temperature reconstructions for all 500 ensemble members of a given method.  
#	Format: Matrix of size 6 (latitude bands) x 121 (ages) x 500 (ensemble members).  
#	Latitude bands are 30-degree bands, ordered from north to south. 
#	latband_ranges: Descriptions of the latitude ranges encompassed by each latitude band. 
#	latband_weights: The proportion of area encompassed by each latitude band, assuming Earth is a perfect sphere.  
#	These are used as weights when computing area-weighted global average of the six latitude bands.
#	
#	temp12k_[method]_globalmean.csv: These are CSV files containing global-mean temperature reconstructions 
#	for all 500 ensemble members of a given method, as in the netCDF file above.  The first column is age, 
#	followed by the global-mean temperature reconstruction for all 500 ensemble members.
#	
#	temp12k_[method]_latband_[region].csv. These files are structured the same way as the "globalmean" files above, 
#	but contain values for latitude-zone means.
#	
#	temp12k_allmethods_percentiles.csv: This is a CSV file containing the multi-method median values as well as 
#	5th and 95th percentiles for the global mean and latitudinal-zone means.  Quantities are calculated over the 
#	multi-method 2500-member ensemble and are not separated by method.  The first column is age and the following 
#	columns are the 5th percentile, median, and 95th percentile for the global mean and each of the six 
#	latitude-zone means. Labels for the columns are given in the first row.
#	
#	ADDITIONAL DETAILS: 
#	Units: Units for ages are years BP, with "present" defined as 1950 CE.  
#	Age values represent the center of 100-year age bins, so a value of 100 yr BP represents the interval 
#	150 to 50 yr BP or 1800-1900 CE.  Units for the temperature reconstructions are anomalies in degrees C. 
#	
#	Temperature reference period: For each of the 2500 ensemble members, the mean over the entire 12000-year 
#	period was removed. Then, for each method, the median of that method's 500 ensemble members was set to 
#	0 degrees C at the time 100 yr BP (1800-1900 CE). This means that the general reference period is 1800-1900 CE, 
#	but uncertainties are calculated at every time step. 
#	
#	The code used to make the calculations can be obtained at: 
#	https://zenodo.org/record/3888590#.XuF2zy2ZODU
#	
#--------------------
# Publication 
#	Authors: Darrell Kaufman, Nicholas McKay, Cody Routson, Michael Erb, Christoph Dätwyler, Philipp Sommer, Oliver Heiri, Basil Davis 
#	Published_Date_or_Year: 2020-06-30
#	Published_Title: Holocene global mean surface temperature, a multi-method reconstruction approach
#	Journal_Name: Scientific Data
#	Volume: 
#	Edition: 
#	Issue: 
#	Pages:
#	Report_Number: 
#	DOI: 10.1038/s41597-020-0530-7
#	Online_Resource:
#	Full_Citation: 
#	Abstract: An extensive new multi-proxy database of paleo-temperature time series (Temperature 12k) enables more robust analysis of global mean surface temperature (GMST) and associated uncertainties than was previously available. We applied five different statistical methods to reconstruct the GMST of the past 12,000 years (Holocene). Each method used different approaches to averaging the globally distributed time series and to characterizing various sources of uncertainty, including proxy temperature, chronology and methodological choices. The results were aggregated to generate a multi-method ensemble of plausible GMST and latitudinal-zone temperature reconstructions with a realistic range of uncertainties. The warmest 200-year-long interval took place around 6500 years ago when GMST was 0.7C (0.4, 1.8) warmer than the 19th Century (median, 5th, 95th percentiles). Following the global thermal maximum, GMST cooled at an average rate -0.08C per 1000 years (-0.24, -0.05). The multi-method ensembles and the code used to generate them highlight the utility of the Temperature 12k database, and they are now available for future use by studies aimed at understanding Holocene evolution of the Earth system.
#------------------
# Funding_Agency 
#	Funding_Agency_Name: US National Science Foundation
#	Grant: AGS-1602105, AGS-1903548
#------------------
# Funding_Agency 
#	Funding_Agency_Name: Swiss National Science Foundation 
#	Grant: SNF 200021-165494
#------------------
# Site_Information 
#	Site_Name: Global
#	Location: Geographic Region>Global
#	Country:  
#	Northernmost_Latitude: 90 
# 	Southernmost_Latitude: -90
# 	Easternmost_Longitude: 180
# 	Westernmost_Longitude: -180
# 	Elevation:  
#------------------
# Data_Collection   
#	Collection_Name: Kaufman2020GMST
#	Earliest_Year: 12000
#	Most_Recent_Year: 0 
#	Time_Unit: Cal. Year BP
#	Core_Length: 
#	Notes: 
#------------------
# Chronology_Information
#	Chronology: 
# 	
#----------------
# Variables 
#
# Data variables follow are preceded by "##" in columns one and two.
# Data line variables format:  one per line, shortname-tab-variable components (what, material, error, units, seasonality, data type,detail, method, C or N for Character or Numeric data, free text) 
#
## age_BP	age, , , calendar year before present, , climate reconstructions, , ,N, 100-year bins
## temp	air temperature, multiple proxies, , degree Celsius, annual, climate reconstructions, anomalized,,N, median across ensemble members of all methods; anomalized relative to 1800-1900CE 
## temp5 air temperature, multiple proxies, percentile, degree Celsius, annual, climate reconstructions, anomalized,,N, 5th percentile across ensemble members of all methods; anomalized relative to 1800-1900CE
## temp95 air temperature, multiple proxies, percentile, degree Celsius, annual, climate reconstructions, anomalized,,N, 95th percentile across ensemble members of all methods; anomalized relative to 1800-1900CE
## tempSCC air temperature, multiple proxies, ensemble member, degree Celsius, annual, climate reconstructions, anomalized,,N, Standard calibrated composite method; anomalized relative to 1800-1900CE 
## tempDCC air temperature, multiple proxies, ensemble member, degree Celsius, annual, climate reconstructions, anomalized,,N, Dynamic calibrated composite method; anomalized relative to 1800-1900CE
## tempGAM air temperature, multiple proxies, ensemble member, degree Celsius, annual, climate reconstructions, anomalized,,N, Generalized additive model method; anomalized relative to 1800-1900CE
## tempCPS air temperature, multiple proxies, ensemble member, degree Celsius, annual, climate reconstructions, anomalized,,N, Composite plus scale method; anomalized relative to 1800-1900CE
## tempPAI air temperature, multiple proxies, ensemble member, degree Celsius, annual, climate reconstructions, anomalized,,N, Pairwise comparison (PaiCo) method; anomalized relative to 1800-1900CE
#
#----------------
# Data:
# Data lines follow (have no #)
# Data line format - tab-delimited text, variable short name as header
# Missing Values:  
#
age_BP	temp
