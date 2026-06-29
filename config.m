
%% Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
%  Configuration file to be read by main scripts, for configuring filepaths
%  to input/output data, centerline search parameters and other user 
%  specifiable parameters. Update parameters as required! 
% 
%  (c) Dylan Kreynen
%  University of Oslo
%  2024-2026
% 
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% user specifiable variables
%  (update as required)

% ice shelf DEM (should be GeoTIFF)
path_to_DEM = './input/venable.tif'; % to a single .tif for map_multi_channel.m, a directory of .tifs for map_channel_timeseries.m
DEM_nodata = -9999;         % DEM no data value
window_DEM = 0;             % window size for DEM smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)
                            % only works with map_multi_channel.m variations (TO DO: implement in map_channel_timeseries.m)

% output behaviour: 
results_dir = './output/';
proj_subdir = 'venable/'; 
fig_subdir = 'fig/'; 
shp_subdir = 'shp/'; 
file_prefix = 'default_';   % (optional)
% output path will be constructed as follows: 
% results_dir/proj_subdir/fig_subdir/file_prefix_....ext

save_figs = 0;              % print figures to disk Y/N
figs_filetype = '-dpng';    % for use with "print()"
figs_resolution = '-r500';  % for use with "print()"
ext_figs = 1;               % plot (and print) extended figures Y/N
save_shps = 0;              % save output as shapefiles Y/N

% select method to specify channel start/end points
start_end_method = 2;
% 1 = click on start/end points
% 2 = read from shapefile
% 3 = manually enter in script
path_to_start_end_shp = './input/venable_start_end_4.shp'; 
% ^ only needed when start_end_method is set to "2" (read from shapefile)
% important: only works if shapefile has same map projection as DEM! 
shelf_filter = ''; % no filtering: shelf_filter = ''; 
% ^ optional: filter shapefile on "Shelf" field

% centerline search parameters
search_step = 1500;               % distance to step away from last known centerline point to construct search profile [m]
search_angle = 85;                % angle of view within to look for next centerline point [deg]
max_gradient = 2;                 % if new centerline point's elevation exceeds max. gradient, pick next best point instead [%]
window_cent = 20;                 % window size for search profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)
max_recursions = 1;               % keep trying with slightly different search parameters in case channel end point is not found [-] (set to 0 or 1 for no recursion)
% length of search segment = (search_angle/360)*2*pi*search_step

% channel cross sectional profile parameters
prof_length   = 6000;             % length of cross sectional profiles [m]
prof_interval = 500;                % spacing between profiles along centerline [m]
                                  % (0 = one profile per centerline search segment)
% channel edge parameters
edge_method = "KneePoint";        % method to use to identify channel edges ("SlopeThreshold", "KneePoint" or "NearPeaks")
knee_method = "LinearRegression"; % kneepoint algorithm to use ("LinearRegression" or "Kneedle", used with edge methods "KneePoint" and "NearPeaks")
min_width = 500;                  % minimum channel width [m] (set to 0 for no minimum width)
max_width = 0;                    % maximum channel width [m] (set to 0 for maximum width = prof_length)
sg_window = 1000;                 % window size for profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing, Savitzky-Golay filter)
m_window = 0;                     % window size for edge smoothing [-] (no. of profile edges, set to 0 for no smoothing, median filter)

slope_thr = 0.00;                 % slope threshold for identifying edge [deg] (only used when edge_method = "SlopeThreshold")
peak_prom = 0.5;                  % MinPeakProm for findpeaks() [m], minimum prominence for channel edge (only used when edge_method = "NearPeaks")
keep_pks = 0;                     % prevent peaks from being adjusted by along-channel edge smoothing (0 or 1, only used when edge_method = "NearPeaks")

% outlier filtering of edges (based on z-score)
z_thr_elev = [1.5 1.5];           % z-score outlier threshold, elevation of edge (set to 0 to skip outlier identification)
                                  % optional: [left_thr right_thr] to use different thresholds for the left and right edge
z_thr_idx = [1.5 1.5];            % z-score outlier threshold, profile index of edge (set to 0 to skip outlier identification)
                                  % optional: [left_thr right_thr] to use different thresholds for the left and right edge
edge_subst_window = 5;            % window size for outlier substitution (moving median filter, set to 0 to leave out outliers altogether, default: 5)
% note: set m_window = 0 (no edge smoothing) when filtering outliers

% some controls for the overview figure, related to edge outlier filtering
plot_prof_transects = 1;          % plot profile transects on overview figure Y/N
plot_edge_gaps = 1;               % plot across the gaps left after edge outlier filtering Y/N


%% manage directories
%  (no need to update)

% add "functions" directory to search path
addpath("./functions")

% create output directories, if necessary
if save_figs % figures
    fig_dir = append(results_dir, proj_subdir, fig_subdir); 
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir)
    end
end

if save_shps % shapefiles
    shp_dir = append(results_dir, proj_subdir, shp_subdir); 
    if ~exist(shp_dir, 'dir')
        mkdir(shp_dir)
    end
end