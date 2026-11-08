
%% Sub-Ice: DEM-based semi-automated mapping of ice shelf basal channels
%  Configuration file to be read by the main scripts (map_multi_channel.m,
%  map_channel_timeseries.m). Update settings as required!
%
%  Settings are grouped by how often you will need to change them:
%   1. input and output: set for every run
%   2. centerline search: tune for your DEM and channels
%   3. profiles and channel edges: defaults as used in the Sub-Ice
%      publication, adjust to the size of your channels if needed
%   4. other settings: rarely changed
%  See the Sub-Ice publication for details. Settings left out of this file
%  are set to their default value (with a warning), see
%  functions/check_config.m.
%
%  (c) Dylan Kreynen
%  University of Oslo
%  2024-2026
%
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% 1. INPUT AND OUTPUT
%  (set for every run)

% ice shelf DEM (GeoTIFF, projected coordinates in meters)
path_to_DEM = './input/venable.tif';    % a single .tif for map_multi_channel.m, a directory of .tifs for map_channel_timeseries.m
                                        % (time series: DEM filenames should end in the acquisition date, yyyymmdd, e.g. venable_20200115.tif)
DEM_nodata = -9999;                     % DEM no data value (read from the GeoTIFF where given; this value is used as fallback)

% channel start and end points
start_end_method = 2;
% 1 = click on start/end points (on a map of the DEM)
% 2 = read from multipoint shapefile (path_to_start_end_shp)
% 3 = give coordinates below (start_end_points)
% 4 = load existing centerlines from line shapefile (path_to_centerline_shp, skips the centerline search)
path_to_start_end_shp = './input/venable_start_end_4.shp';
% ^ start_end_method = 2: points ordered start, end, start, end, ... (one
%   pair per channel; time series: one pair per DEM, in chronological order)
%   important: only works if shapefile has same map projection as DEM!
start_end_points = [];
% ^ start_end_method = 3: [x_start y_start x_end y_end] in map coordinates
%   (same projection as DEM), one row per channel (time series: per DEM)
path_to_centerline_shp = './output/venable/shp/default_all_centerlines.shp';
% ^ start_end_method = 4: line shapefile, e.g. a centerline selected from a
%   previous (ensemble) run; lines should run from channel start to end point
%   (downstream); channel labels are read from the first attribute starting
%   with "channel" (if any), other attributes are carried over to the output
%   important: only works if shapefile has same map projection as DEM!

% output
% (output path: results_dir/proj_subdir/fig_subdir/file_prefix....ext)
results_dir = './output/';
proj_subdir = 'venable/';
file_prefix = 'default_';               % (optional)
save_figs = 0;                          % save figures Y/N
save_shps = 0;                          % save centerlines, channel edges and profiles as shapefiles Y/N
save_tables = 0;                        % save profiles and channel metrics as tables (.csv) Y/N
ext_figs = 1;                           % extended figures: cross-sectional profiles and metrics along each channel Y/N
                                        % (shown when not saved, saved without opening windows when save_figs = 1)


%% 2. CENTERLINE SEARCH
%  tune for your DEM and channels: the search depends on DEM resolution,
%  channel size and surface character; if centerlines do not follow the
%  channel or do not reach the end point, try other values, or ensemble mode

% single values: map one centerline per channel (and continue with profiles
% and edges); sets of values for any of search_step, search_angle,
% max_gradient, window_cent and window_DEM: ensemble mode, map one
% centerline per channel for every combination of values (full-factorial
% ensemble, e.g. search_angle = 30:10:120), then stop - select a centerline
% per channel from the ensemble output and continue with start_end_method = 4
window_DEM = 0;                         % window size for DEM smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)
search_step = 1000;                     % distance to step away from last known centerline point to construct search profile [m]
search_angle = 85;                      % angle of view within to look for next centerline point [deg]
max_gradient = 2;                       % if new centerline point's elevation exceeds max. gradient, pick next best point instead [%]
window_cent = 20;                       % window size for search profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)
max_recursions = 10;                    % max. no. of centerline search attempts in case channel end point is not found [-] (1 = no retries, 2 = retry from end point only,
                                        % > 2 = then keep retrying with randomly perturbed search parameters, see find_centerline())
% length of search segment = (search_angle/360)*2*pi*search_step


%% 3. PROFILES AND CHANNEL EDGES
%  defaults as used in the Sub-Ice publication; adjust prof_length,
%  min_width and max_width to the size of your channels if needed

% channel cross sectional profile parameters
prof_length   = 8000;                   % length of cross sectional profiles [m]
prof_interval = 250;                    % spacing between profiles along centerline [m] (0 = one profile per centerline search segment)

% channel edge parameters
edge_method = "KneePoint";              % method to use to identify channel edges ("SlopeThreshold", "KneePoint" or "NearPeaks")
knee_method = "Kneedle";                % kneepoint algorithm to use ("LinearRegression" or "Kneedle", used with edge methods "KneePoint" and "NearPeaks")
min_width = 1000;                       % minimum channel width [m] (set to 0 for no minimum width; or [left right] min. distance from centerline per side)
max_width = 8000;                       % maximum channel width [m] (set to 0 for maximum width = prof_length; or [left right] max. distance from centerline per side)
sg_window = 1000;                       % window size for profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing, Savitzky-Golay filter)
m_window = 0;                           % window size for edge smoothing [-] (no. of profile edges, set to 0 for no smoothing, median filter)
slope_thr = 0.00;                       % slope threshold for identifying edge [deg] (only used when edge_method = "SlopeThreshold")
peak_prom = 0.25;                       % MinPeakProm for findpeaks() [m], minimum prominence for channel edge (only used when edge_method = "NearPeaks")
keep_pks = 0;                           % prevent peaks from being adjusted by along-channel edge smoothing (0 or 1, only used when edge_method = "NearPeaks")

% outlier filtering of edges (based on z-score)
% note: use either edge smoothing (m_window) or outlier filtering, not both
z_thr_elev = [1.75 1.75];               % z-score outlier threshold, elevation of edge (set to 0 to skip outlier identification)
                                        % optional: [left_thr right_thr] to use different thresholds for the left and right edge
z_thr_idx = [1.75 1.75];                % z-score outlier threshold, profile index of edge (set to 0 to skip outlier identification)
                                        % optional: [left_thr right_thr] to use different thresholds for the left and right edge
edge_subst_window = 5;                  % window size for outlier substitution (moving median filter, set to 0 to leave out outliers altogether, default: 5)

% outlier filtering of trough depth (based on z-score)
z_thr_trough_elev = 0;                  % z-score outlier threshold, elevation of trough (set to 0 to skip outlier identification)
trough_subst_window = 3;                % window for outlier substitution (moving mean filter, set to 0 to leave out outliers altogether, default: 3)


%% 4. OTHER SETTINGS
%  rarely changed

% output details
fig_subdir = 'fig/';
shp_subdir = 'shp/';
tab_subdir = 'tab/';
fig_format = 'png';                     % figure file format (e.g. 'png', 'pdf', 'svg')
fig_dpi = 500;                          % figure resolution [dpi]

% optional channel selection
shelf_filter = '';                      % filter start/end point shapefile on "Shelf" attribute (no filtering: '')
excl_channel_labels = [];               % channel labels to skip, e.g. ["channel-3", "channel-7"] (no exclusion: [])

% map figures
map_hillshade = 1;                      % add a hillshade to the DEM on map figures Y/N
map_clims = [];                         % [min max] elevation for the DEM colors on map figures [m] ([] = automatic, 2nd-98th percentile within the map)
plot_prof_transects = 1;                % plot profile transects on overview figure Y/N
plot_edge_gaps = 1;                     % plot across the gaps left after edge outlier filtering Y/N
