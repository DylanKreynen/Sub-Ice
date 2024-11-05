## Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels

Mapping algorithm that returns ice shelf basal channels' centerline, outlines
and cross sectional profiles; based on surface expressions in a DEM. 
For now, please refer to comments in the code for details! 


 #### Input: 
  - digital elevation model (GeoTIFF)
  - channels' start and end points (through GUI or multipoint shapefile)

 #### Output: 
  - channel centerlines and outlines (figures and georeferenced shapefiles)
  - cross sectional elevation/depth profiles (figures)
  - channel depth and width along centerlines (figures)


#### Getting started: 
1. Clone this GitHub repository to create a local copy of the code on your computer (or download and unzip). 
2. Update the user specifiable parameters in *config.m* as required. 
3. Run *map_multi_channel.m* using a recent version of Matlab (code tested in R2024a). 
That's it! All data required to map your first basal channels is already included. 


 #### Main scripts in this distribution: 
  - *map_multi_channel.m*:
    map one or more channels given a single DEM
  - *map_channel_timeseries.m*:
    map a single channel over different DEMs
  - *config.m*
    configuration file with user specifiable variables

 #### Directories in this distribution: 
  - *./functions*:
    collection of custom Matlab functions required to run main scripts
  - *./input*:
    sample [REMA elevation data](https://www.pgc.umn.edu/data/rema/) and shapefiles with channel start/end points
  - *./output*:
    sample output (will be overwritten when running main scripts)


\
**Sub-Ice**: **Su**rface-**b**ased **I**ce Shelf **C**hannel **E**xtraction  
version as of 2024-07-11 (tested using R2024a)

&copy; [Dylan Kreynen](https://www.mn.uio.no/geo/english/people/aca/geohyd/kreynen/)  
University of Oslo  
2024

originally a project at the International Summer School in Glaciology\
project team members: Marcelo Santis & Dylan Kreynen\
advisor: Karen Alley (University of Manitoba)\
McCarthy (AK), June 2024

![Sub-Ice example output.](output/fig_readme.png?raw=true)