
%% Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
%  Configuration file to be read by main scripts, for configuring filepaths
%  to input/output data, centerline search parameters and other user 
%  specifiable parameters. Update parameters as required! 
% 
%  (c) Dylan Kreynen
%  University of Oslo
%  June - August 2024
% 
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% user specifiable variables
%  (update as required)

% path to DEM (should be GeoTIFF)
path_to_DEM = '..\REMA\getz_200m.tif'; 

% output behaviour: 
results_dir = '.\output\';
proj_subdir = 'getz_200m\1000_20_10\'; 
fig_subdir = 'fig\'; 
shp_subdir = 'shp\'; 
file_prefix = 'default_'; % (optional)
% output path will be constructed as follows: 
% results_dir\proj_subdir\fig_subdir\file_prefix_....ext

save_figs = 1;              % print figures to disk Y/N
figs_filetype = '-dpng';    % for use with "print()"
figs_resolution = '-r500';  % for use with "print()"
ext_figs = 0;               % plot (and print) extended figures Y/N
save_shps = 0;              % save output as shapefiles Y/N

% select method to specify channel start/end points
start_end_method = 2;
% 1 = click on start/end points
% 2 = read from shapefile
% 3 = manually enter in script
path_to_start_end_shp = '.\input\Getz1_2\Main_channels.shp'; 
%path_to_start_end_shp = '.\input\Getz1_2\Short_channels.shp'; 
%path_to_start_end_shp = '.\input\Getz1_2\Tributaries.shp'; 
% ^ only needed when start_end_method is set to "2" (read from shapefile)
% important: only works if shapefile has same map projection as DEM! 

% centerline search parameters
search_step = 1000;               % distance to step away from last known centerline point to construct search profile [m]
no_cent_samp_pts = 20;            % number of sampling points on search profile [-]
cent_samp_step = 100;              % distance between sampling points on search profile [m]
max_no_cent_pts = 70;             % when to stop looking for centerline end point [-]
crack_thr = 10;                    % if new centerline point's depth w.r.t. last known point is greater than threshold, pick next best point instead [m]
window_cent = 0;                  % window size for search profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)

% channel cross sectional profile parameters
prof_samp_step = 50;              % sdistance between sampling points on profile [m]
no_prof_samp_pts = 100;           % number of sampling points on profile [-]
% note: profile length ~ no_sampling_points*prof_samp_step

% channel edge parameters
slope_thr = 0.25;                 % slope threshold for identifying edge [deg]
window_edge = 500;                % window size for profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)

% add "functions" directory to search path
addpath("./functions")