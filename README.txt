Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
ver 2024-07-08 (tested in Matlab 2024a, requires the Mapping Toolbox)

Algorithm that returns ice shelf basal channels' centerline, outlines
and cross sectional profiles; based on surface expressions in a DEM. 
See comments in code for details. 

 Input: 
  - DEM (GeoTIFF)
  - channel start and end points (through GUI or read from shapefile)
 
 Output: 
  - channel centerlines and outlines (figures and georeferenced shapefiles)
  - cross sectional elevation/depth profiles (figures)
  - channel depth and width along centerlines (figures)
 
 
 Main scripts in this distribution: 
  - map_multi_channel.m
    map one or more channels given a single DEM
  - map_channel_timeseries.m
    map single channel over different DEMs
  
 Directories in this distribution: 
  - ./functions
    collection of custom Matlab functions required to run main scripts
  - ./input
    sample REMA elevation data and shapefiles with channel start/end pts
  - ./output
    sample output (will be overwritten when running main scripts)


Sub-Ice: SUrface-based Ice shelf Channel Extraction scheme
 
(c) Dylan Kreynen
University of Oslo
June - July 2024
 
originally a project at the Int. Summer School in Glaciology
project team members: Marcelo Santis & Dylan Kreynen
advisor: Karen Alley (University of Manitoba)
McCarthy (AK), June 2024