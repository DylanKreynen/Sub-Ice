%% Sub-Ice: check configuration
%  Script run by the main scripts right after the configuration file:
%   - settings missing from the configuration file are set to their
%     default value (listed below), with a warning
%   - settings are checked for obvious mistakes: clear error messages
%     where Sub-Ice cannot run, warnings where it can
%   - output directories are created (if output is switched on)
%  No need to update; change settings in the configuration file instead.
%
%  (c) Dylan Kreynen
%  University of Oslo
%  2026


%% default values
%  (as in config.m; used for settings missing from the configuration file)

cfg_defaults = {
    % input and output
    'DEM_nodata',           -9999
    'window_DEM',           0
    'results_dir',          './output/'
    'proj_subdir',          ''
    'fig_subdir',           'fig/'
    'shp_subdir',           'shp/'
    'tab_subdir',           'tab/'
    'file_prefix',          ''
    'save_figs',            0
    'fig_format',           'png'
    'fig_dpi',              500
    'ext_figs',             0
    'save_shps',            0
    'save_tables',          0
    % channel start and end points
    'start_end_method',     2
    'path_to_start_end_shp', ''
    'start_end_points',     []
    'path_to_centerline_shp', ''
    'shelf_filter',         ''
    'excl_channel_labels',  []
    % centerline search
    'search_step',          1500
    'search_angle',         85
    'max_gradient',         2
    'window_cent',          20
    'max_recursions',       1
    % profiles
    'prof_length',          8000
    'prof_interval',        0
    % edges
    'edge_method',          "KneePoint"
    'knee_method',          "Kneedle"
    'min_width',            1000
    'max_width',            8000
    'sg_window',            500
    'm_window',             0
    'slope_thr',            0.00
    'peak_prom',            0.25
    'keep_pks',             0
    % outliers
    'z_thr_elev',           [0 0]
    'z_thr_idx',            [1.75 1.75]
    'edge_subst_window',    5
    'z_thr_trough_elev',    0
    'trough_subst_window',  3
    % maps
    'map_hillshade',        1
    'map_clims',            []
    'plot_prof_transects',  1
    'plot_edge_gaps',       1
    };

% (older configuration files: figure format and resolution given for print())
if exist('figs_filetype', 'var') && ~exist('fig_format', 'var')
    fig_format = erase(string(figs_filetype), "-d");
end
if exist('figs_resolution', 'var') && ~exist('fig_dpi', 'var')
    fig_dpi = str2double(erase(string(figs_resolution), "-r"));
end

% set missing settings to their default value
cfg_missing = strings(0);
for cfg_k = 1:size(cfg_defaults, 1)
    if ~exist(cfg_defaults{cfg_k, 1}, 'var')
        eval(append(cfg_defaults{cfg_k, 1}, ' = cfg_defaults{cfg_k, 2};'))
        cfg_missing(end+1) = cfg_defaults{cfg_k, 1}; %#ok<SAGROW>
    end
end
% (no warning for settings not used with the chosen start/end method, nor
% for the optional channel selection)
cfg_silent = ["shelf_filter", "excl_channel_labels"];
if start_end_method ~= 2, cfg_silent(end+1) = "path_to_start_end_shp"; end
if start_end_method ~= 3, cfg_silent(end+1) = "start_end_points"; end
if start_end_method ~= 4, cfg_silent(end+1) = "path_to_centerline_shp"; end
cfg_missing = setdiff(cfg_missing, cfg_silent, 'stable');
if ~isempty(cfg_missing)
    disp(append("Warning: not set in the configuration file, using default values for: ", strjoin(cfg_missing, ", "), ". "))
end


%% checks
%  errors where Sub-Ice cannot run, warnings where it can

% DEM (single file or, for a time series, a directory)
if ~exist('path_to_DEM', 'var')
    error("No DEM given: set path_to_DEM in the configuration file. ")
end
if exist('subice_mode', 'var') && subice_mode == "timeseries"
    if ~isfolder(path_to_DEM)
        error(append("path_to_DEM should be a directory with DEMs for map_channel_timeseries.m (now: ", path_to_DEM, "). "))
    end
elseif isfolder(path_to_DEM)
    error(append("path_to_DEM is a directory (", path_to_DEM, "): give a single DEM for map_multi_channel.m, or use map_channel_timeseries.m for a time series of DEMs. "))
elseif ~isfile(path_to_DEM)
    error(append("DEM not found: ", path_to_DEM, " (check path_to_DEM in the configuration file). "))
end

% start and end points
if ~ismember(start_end_method, 1:4)
    error(append("start_end_method should be 1 (click), 2 (shapefile), 3 (coordinates in configuration file) or 4 (centerlines from shapefile), not ", string(start_end_method), ". "))
end
if start_end_method == 2 && ~isfile(path_to_start_end_shp)
    error(append("Start/end point shapefile not found: '", path_to_start_end_shp, "' (check path_to_start_end_shp in the configuration file). "))
end
if start_end_method == 3 && (isempty(start_end_points) || size(start_end_points, 2) ~= 4)
    error("start_end_method = 3: give the start and end points in the configuration file as start_end_points = [x_start y_start x_end y_end], one row per channel (map coordinates). ")
end
if start_end_method == 4 && ~isfile(path_to_centerline_shp)
    error(append("Centerline shapefile not found: '", path_to_centerline_shp, "' (check path_to_centerline_shp in the configuration file). "))
end

% edge detection method (not case sensitive, spaces ignored)
cfg_edge_methods = ["SlopeThreshold", "KneePoint", "NearPeaks"];
cfg_k = find(strcmpi(regexprep(string(edge_method), '\s', ''), cfg_edge_methods));
if isempty(cfg_k)
    error(append("edge_method should be ", strjoin(cfg_edge_methods, ", "), ", not '", string(edge_method), "'. "))
end
edge_method = cfg_edge_methods(cfg_k);
cfg_knee_methods = ["LinearRegression", "Kneedle"];
cfg_k = find(strcmpi(regexprep(string(knee_method), '\s', ''), cfg_knee_methods));
if isempty(cfg_k)
    error(append("knee_method should be ", strjoin(cfg_knee_methods, " or "), ", not '", string(knee_method), "'. "))
end
knee_method = cfg_knee_methods(cfg_k);

% channel width limits
if all(max_width > 0) && any(sum(min_width) >= sum(max_width))
    error("min_width should be smaller than max_width (check the configuration file). ")
end
if all(max_width > 0) && sum(max_width) > prof_length
    disp(append("Warning: max_width (", string(sum(max_width)), " m) exceeds prof_length (", string(prof_length), " m): channel width is limited by the profile length. "))
end

% edge smoothing and outlier filtering: use either, not both
if m_window > 0 && (any(z_thr_idx > 0) || any(z_thr_elev > 0))
    disp("Warning: both edge smoothing (m_window) and edge outlier filtering (z_thr_idx, z_thr_elev) are switched on. Smoothing is applied first; consider using either, not both. ")
end

% figure output (as used with print())
figs_filetype = append('-d', char(fig_format));
figs_resolution = append('-r', num2str(fig_dpi));


%% output directories
%  results_dir/proj_subdir/fig_subdir etc. (created if output is switched on)

if save_figs
    fig_dir = append(results_dir, proj_subdir, fig_subdir);
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir), end
end
if save_shps
    shp_dir = append(results_dir, proj_subdir, shp_subdir);
    if ~exist(shp_dir, 'dir'), mkdir(shp_dir), end
end
if save_tables
    tab_dir = append(results_dir, proj_subdir, tab_subdir);
    if ~exist(tab_dir, 'dir'), mkdir(tab_dir), end
end

clear cfg_defaults cfg_missing cfg_silent cfg_k cfg_edge_methods cfg_knee_methods
