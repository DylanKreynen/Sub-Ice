clear all
close all
clc

%% Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
%  Algorithm that returns ice shelf basal channels' centerline, outlines
%  and cross sectional profiles; based on surface expressions in a DEM.
%  This script maps a single channel over a time series of DEMs, given
%  user identified channel start and end points on each DEM, or existing
%  channel centerlines (one per DEM). (See "map_multi_channel.m" for
%  mapping one or more channels given a single DEM.)
%
%  The DEMs are read from a directory (path_to_DEM in config.m); the
%  acquisition date of each DEM is read from its filename, which should
%  end in the date as yyyymmdd (e.g. "venable_20200115.tif"). DEMs are
%  processed in chronological order. All other settings are as for
%  map_multi_channel.m (single centerline per DEM: no ensemble mode).
%
%  Input:
%   - DEMs (GeoTIFFs, in one directory, date in filename)
%   - channel start and end points, one pair per DEM in chronological
%     order (through GUI, read from shapefile or entered manually), or
%     channel centerlines, one per DEM in chronological order (shapefile)
%
%  Output:
%   - channel centerlines and outlines per date (figures and georeferenced shapefiles)
%   - cross sectional elevation/depth profiles per date (figures and tables)
%   - channel depth and width along centerlines per date (figures and tables)
%   - summary of channel metrics per date (table)
%
%  (c) Dylan Kreynen
%  University of Oslo
%  2024-2026
%
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% run configuration file - set user specifiable parameters
%  contains filepaths to input/output data, sets method to select channel
%  start/end points, configures centerline search parameters, output
%  behaviour etc. > update configuration file as required!
%  (default: config.m; when run from run_batch.m: the configuration file
%  passed on by run_batch.m)

config_fn = 'config.m';
if isappdata(groot, 'subice_config_fn')
    config_fn = getappdata(groot, 'subice_config_fn');
    rmappdata(groot, 'subice_config_fn')    % (applies to this run only)
end
disp(append("Running Sub-Ice (time series) with configuration file ", config_fn, ". "))
% (evaluated in the current folder: relative paths in the configuration
% file are relative to the Sub-Ice folder, wherever the file is stored)
eval(fileread(config_fn))

% check configuration: default values for missing settings, clear messages
% for mistakes, create output directories (see functions/check_config.m)
addpath("./functions")
subice_mode = "timeseries";
check_config

% no ensemble mode for time series
if any(cellfun(@numel, {window_DEM, search_step, search_angle, max_gradient, window_cent}) > 1)
    error("Ensemble mode (sets of centerline search parameters in config.m) is not available for time series: give single values. ")
end

% figure style (colors, font sizes; see fig_style())
sty = fig_style();


%% list DEMs/GeoTIFFs
%  and corresponding dates (from filename, ending in yyyymmdd)

if endsWith(path_to_DEM, ".tif", 'IgnoreCase', true) || ~isfolder(path_to_DEM)
    error("path_to_DEM should point to a directory with DEMs for map_channel_timeseries.m (check in config.m).")
end

DEM_files = dir(fullfile(path_to_DEM, '*.tif'));
DEM_files = {DEM_files.name}';
no_DEMs = length(DEM_files);
if no_DEMs == 0
    error(append("No DEMs (.tif) found in ", path_to_DEM, ". "))
end

dates = NaT(no_DEMs, 1);
for t = 1:no_DEMs
    [~, name] = fileparts(DEM_files{t});
    yyyymmdd = regexp(name, '\d{8}$', 'match', 'once');
    if isempty(yyyymmdd)
        error(append("Could not read date from DEM filename ", DEM_files{t}, " (filename should end in yyyymmdd). "))
    end
    dates(t) = datetime(yyyymmdd, 'InputFormat', 'yyyyMMdd');
end

% sort chronologically
[dates, idx] = sort(dates);
DEM_files = DEM_files(idx);
date_label = string(dates, 'yyyy-MM-dd');
disp(append("Found ", string(no_DEMs), " DEMs, from ", date_label(1), " to ", date_label(end), ". "))


%% specify channel centerline start/end points
%  have user click on channel start and end points through GUI (on each
%  DEM), read start and end points from shapefile or enter img coordinates
%  manually; or load existing channel centerlines from shapefile
%  (start/end points and centerlines are kept in map coordinates, as the
%  DEMs may differ in extent)

attr_startend = table('Size', [no_DEMs 0], 'VariableTypes', {});

if start_end_method == 1
    % click on start/end points, on each DEM
    P_start_w = NaN(no_DEMs, 2);
    P_end_w = NaN(no_DEMs, 2);
    for t = 1:no_DEMs
        [DEM, R] = read_dem(fullfile(path_to_DEM, DEM_files{t}), DEM_nodata);
        DEM = double(DEM);
        figure(1)
        clf
        plot_map_base(gca, DEM, R, 'clims', map_clims, 'hillshade', map_hillshade);
        title(append(date_label(t), ": left click on channel start and end point, respectively. "))
        disp(append(date_label(t), ": left click on channel start and end point, respectively. "))
        [xs, ys] = ginput(1);
        scatter(xs, ys, 'm', 'filled')
        [xe, ye] = ginput(1);
        scatter(xe, ye, 'm', 'filled')
        [P_start_w(t,1), P_start_w(t,2)] = intrinsicToWorld(R, xs, ys);
        [P_end_w(t,1), P_end_w(t,2)] = intrinsicToWorld(R, xe, ye);
        pause(0.5)
    end
    close(figure(1))

elseif start_end_method == 2
    % read start/end points from shapefile: pairs of start and end point,
    % one pair per DEM, in chronological order
    % only works if shapefile has same map projection as the DEMs!
    S = shaperead(path_to_start_end_shp);
    % filter on "shelf" field
    if ~isempty(shelf_filter)
        S = S(strcmp({S.Shelf}, shelf_filter));
    end
    x_startend = vertcat(S.X);
    y_startend = vertcat(S.Y);
    if rem(length(x_startend), 2) ~= 0
        error('Found odd number of start and end points in shapefile. Please provide an even number of points, ordered channel start, end respectively. ');
    end
    even_idx = rem(1:length(x_startend), 2) == 0;
    P_start_w = [x_startend(~even_idx), y_startend(~even_idx)];
    P_end_w = [x_startend(even_idx), y_startend(even_idx)];

    % attributes of start points, carried over to the centerline shapefile
    if size(P_start_w, 1) == no_DEMs
        attr_startend = struct2table(S(~even_idx), 'AsArray', true);
        attr_startend = removevars(attr_startend, intersect(attr_startend.Properties.VariableNames, {'Geometry', 'BoundingBox', 'X', 'Y'}));
    end

elseif start_end_method == 3
    % start/end points given in the configuration file (start_end_points,
    % map coordinates: [x_start y_start x_end y_end], one row per DEM, in
    % chronological order)
    P_start_w = start_end_points(:, 1:2);
    P_end_w = start_end_points(:, 3:4);

elseif start_end_method == 4
    % load existing centerlines from line shapefile, one per DEM, in
    % chronological order (no centerline search)
    % only works if shapefile has same map projection as the DEMs!
    S = shaperead(path_to_centerline_shp);
    if length(S) ~= no_DEMs
        error(append("Found ", string(length(S)), " centerlines in shapefile, but ", string(no_DEMs), " DEMs. Provide one centerline per DEM, in chronological order. "))
    end
    xw_load = cell(no_DEMs, 1);
    yw_load = cell(no_DEMs, 1);
    for t = 1:no_DEMs
        valid = ~isnan(S(t).X) & ~isnan(S(t).Y);
        xw_load{t} = S(t).X(valid)';
        yw_load{t} = S(t).Y(valid)';
    end
    P_start_w = cell2mat(cellfun(@(x, y) [x(1), y(1)], xw_load, yw_load, 'UniformOutput', false));
    P_end_w = cell2mat(cellfun(@(x, y) [x(end), y(end)], xw_load, yw_load, 'UniformOutput', false));
    attr_startend = struct2table(S, 'AsArray', true);
    attr_startend = removevars(attr_startend, intersect(attr_startend.Properties.VariableNames, {'Geometry', 'BoundingBox', 'X', 'Y'}));
    % (date, DEM file and label are written anew, e.g. when loading the
    % output of a previous run; other attributes are carried over)
    vars = attr_startend.Properties.VariableNames;
    attr_startend = removevars(attr_startend, vars(ismember(vars, {'date', 'dem_file'}) | startsWith(vars, "channel", 'IgnoreCase', true)));

else
    error('Specify a valid method to enter channel start/end points (under "user specifiable variables" in config.m). ')
end

% check if we have start/end points for all DEMs
if size(P_start_w, 1) ~= no_DEMs
    error(append("Number of start/end point pairs (", string(size(P_start_w, 1)), ") does not match number of DEMs (", string(no_DEMs), "). "))
end

% check start and end points against each DEM: skip dates with points
% outside the DEM, warn for a different projection (see check_points())
shp_in = "";
if start_end_method == 2, shp_in = path_to_start_end_shp; end
if start_end_method == 4, shp_in = path_to_centerline_shp; end
keep = false(no_DEMs, 1);
for t = 1:no_DEMs
    R_t = georasterinfo(fullfile(path_to_DEM, DEM_files{t})).RasterReference;
    [xs, ys] = worldToIntrinsic(R_t, P_start_w(t,1), P_start_w(t,2));
    [xe, ye] = worldToIntrinsic(R_t, P_end_w(t,1), P_end_w(t,2));
    inside = @(x, y) x >= 0.5 & x <= R_t.RasterSize(2) + 0.5 & y >= 0.5 & y <= R_t.RasterSize(1) + 0.5;
    keep(t) = inside(xs, ys) & inside(xe, ye);
    if t == 1
        check_points(R_t, xs, ys, xe, ye, date_label(t), 'shp', shp_in, 'what', "date");   % (projection check, once)
    end
end
if ~any(keep)
    error("All start and end points lie outside their DEMs. Check that they are in the same projection as the DEMs (map coordinates in meters). ")
elseif any(~keep)
    disp(append("Warning: start or end point outside the DEM for date ", strjoin(date_label(~keep), ", "), ": skipped. "))
    DEM_files = DEM_files(keep);
    dates = dates(keep);
    date_label = date_label(keep);
    P_start_w = P_start_w(keep, :);
    P_end_w = P_end_w(keep, :);
    if height(attr_startend) == numel(keep), attr_startend = attr_startend(keep, :); end
    if start_end_method == 4, xw_load = xw_load(keep); yw_load = yw_load(keep); end
    no_DEMs = sum(keep);
end


%% map channel geometries and extract profiles
%  one DEM at a time; geometries are kept in map coordinates (xw_..., yw_...)

xw_cent = cell(no_DEMs, 1);
yw_cent = cell(no_DEMs, 1);
xw_prof = cell(no_DEMs, 1);
yw_prof = cell(no_DEMs, 1);
edge_coord_w = cell(no_DEMs, 1);
channel_length = cell(no_DEMs, 1);
profiles = cell(no_DEMs, 1);
edge_idx = cell(no_DEMs, 1);
edge_elev = cell(no_DEMs, 1);
valid_edges = cell(no_DEMs, 1);
trough_elev = cell(no_DEMs, 1);
valid_trough = cell(no_DEMs, 1);
T_along = cell(no_DEMs, 1);
cent_info = cell(no_DEMs, 1);
res_t = NaN(no_DEMs, 1);
channel_status = zeros(no_DEMs, 1);

% edge detection method, for figure subtitles
edge_method_txt = string(edge_method);
if edge_method ~= "SlopeThreshold"
    edge_method_txt = append(edge_method_txt, "-", knee_method);   % (knee point method also used by NearPeaks)
end

for t = 1:no_DEMs

    % read DEM from GeoTIFF, replace no data with NaN
    fn_DEM = fullfile(path_to_DEM, DEM_files{t});
    [DEM, R] = read_dem(fn_DEM, DEM_nodata);
    DEM = double(DEM);
    res = R.CellExtentInWorldX;     % resolution of DEM [m]
    res_t(t) = res;

    if start_end_method == 4
        % centerline loaded from shapefile
        [x_cent, y_cent] = worldToIntrinsic(R, xw_load{t}, yw_load{t});
        cent_length = sqrt(diff(x_cent).^2 + diff(y_cent).^2);
        channel_status(t) = 1;
    else
        % smooth DEM for centerline mapping, if requested (cached to disk)
        DEM_smooth = smooth_dem(DEM, res, window_DEM, fn_DEM);

        % find centerline (attention: possibly using smoothed DEM!)
        [x_s, y_s] = worldToIntrinsic(R, P_start_w(t,1), P_start_w(t,2));
        [x_e, y_e] = worldToIntrinsic(R, P_end_w(t,1), P_end_w(t,2));
        [x_cent, y_cent, cent_length, cent_info{t}] = find_centerline([x_s, y_s], [x_e, y_e], DEM_smooth, R, ...
                                                    'search_step',      search_step, ...
                                                    'search_angle',     search_angle, ...
                                                    'max_gradient',     max_gradient, ...
                                                    'window',           window_cent, ...
                                                    'max_recursions',   max_recursions, ...
                                                    'verbose',          false);
        clear DEM_smooth

        % did we reach the channel end point?
        % (if not, continue with the centerline of the last attempt)
        channel_status(t) = cent_info{t}.end_found;
    end
    channel_length{t} = sum(cent_length)*res;   % [m]

    % find cross sectional profiles
    [profiles{t}, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, ...
                                                'prof_length',      prof_length, ...
                                                'prof_interval',    prof_interval);

    % find channel edges/outlines
    [edge_idx{t}, edge_coord, edge_elev{t}, valid_edges{t}] = find_edges(profiles{t}, x_prof, y_prof, res, ...
                                                'edge_method',      edge_method, ...
                                                'knee_method',      knee_method, ...
                                                'min_width',        min_width, ...
                                                'max_width',        max_width, ...
                                                'sg_window',        sg_window, ...
                                                'm_window',         m_window, ...
                                                'slope_thr',        slope_thr, ...
                                                'peak_prom',        peak_prom, ...
                                                'keep_pks',         keep_pks, ...
                                                'z_thr_elev',       z_thr_elev, ...
                                                'z_thr_idx',        z_thr_idx, ...
                                                'edge_subst_window', edge_subst_window, ...
                                                'verbose',          t == 1);     % (messages once, not per date)

    % find channel trough elevation and depth
    [trough_elev{t}, ~, valid_trough{t}] = find_trough(profiles{t}, edge_elev{t}, ...
                                                'z_thr',            z_thr_trough_elev, ...
                                                'subst_window',     trough_subst_window);

    % channel metrics per profile (see channel_metrics())
    T_along{t} = channel_metrics(x_cent, y_cent, x_prof, y_prof, ...
        edge_idx{t}, edge_coord, edge_elev{t}, valid_edges{t}, ...
        trough_elev{t}, valid_trough{t}, R, prof_interval);

    % map of this date
    % (per date figures: not shown when saved to file, else docked; see new_figure())
    f = new_figure([], save_figs);
    plot_overview_map(axes(f), DEM, R, {x_cent}, {y_cent}, {edge_coord}, {x_prof}, {y_prof}, date_label(t), ...
        'clims', map_clims, 'hillshade', map_hillshade, 'transects', plot_prof_transects == 1, ...
        'bridge_gaps', plot_edge_gaps == 1, 'title', append("Mapped channel, ", date_label(t)), ...
        'subtitle', append(DEM_files{t}, ", edge detection: ", edge_method_txt))
    if save_figs
        print(f, append(fig_dir, file_prefix, date_label(t), '_mapped_channel'), figs_filetype, figs_resolution)
        close(f)
    end

    % extended figures of this date
    if ext_figs
        prof_dist_vector = ((1:size(profiles{t}, 1)) - (size(profiles{t}, 1)+1)/2)*res;

        f = new_figure([18 10], save_figs);
        plot_profile_figure(axes(f), prof_dist_vector, profiles{t}, edge_idx{t}, edge_elev{t}, ...
            'title', date_label(t), ...
            'subtitle', append("Cross-sectional profiles (", string(size(profiles{t}, 2)), "), looking upstream"), ...
            'min_width', min_width, 'max_width', max_width)
        if save_figs
            print(f, append(fig_dir, file_prefix, date_label(t), '_full_profiles_elev'), figs_filetype, figs_resolution)
            close(f)
        end

        f = new_figure([18 10], save_figs);
        plot_profile_figure(axes(f), prof_dist_vector, profiles{t}, edge_idx{t}, edge_elev{t}, ...
            'relative', true, 'title', date_label(t), ...
            'subtitle', append("Cross-sectional profiles between the channel edges, looking upstream; mean depth below left edge: ", ...
                string(round(mean(T_along{t}.depth_m, 'omitnan'), 1)), " m"), ...
            'min_width', min_width, 'max_width', max_width)
        if save_figs
            print(f, append(fig_dir, file_prefix, date_label(t), '_lim_profiles_depth'), figs_filetype, figs_resolution)
            close(f)
        end

        f = new_figure([18 13], save_figs);
        plot_along_channel_figure(f, T_along{t}, date_label(t))
        if save_figs
            print(f, append(fig_dir, file_prefix, date_label(t), '_along_channel'), figs_filetype, figs_resolution)
            close(f)
        end
    end

    % keep geometries in map coordinates
    [xw_cent{t}, yw_cent{t}] = intrinsicToWorld(R, x_cent, y_cent);
    [xw_prof{t}, yw_prof{t}] = intrinsicToWorld(R, x_prof, y_prof);
    edge_coord_w{t} = NaN(size(edge_coord));
    [edge_coord_w{t}(:,1), edge_coord_w{t}(:,2)] = intrinsicToWorld(R, edge_coord(:,1), edge_coord(:,2));
    [edge_coord_w{t}(:,3), edge_coord_w{t}(:,4)] = intrinsicToWorld(R, edge_coord(:,3), edge_coord(:,4));

    % progress: one line per date
    msg = append("Mapped ", date_label(t), " (", string(t), "/", string(no_DEMs), ", ", DEM_files{t}, "): centerline ", ...
        string(round(channel_length{t}/1000, 1)), " km, ", string(size(profiles{t}, 2)), " profiles. ");
    if ~channel_status(t)
        msg = append(msg, "Warning: end point not reached, continued with the centerline of the last attempt. ");
    end
    disp(msg)

    pause(0.01) % (force matlab to plot)
end

disp(append("Finished mapping all dates. End point reached for ", string(sum(channel_status)), "/", string(no_DEMs), " dates. "))
if any(~channel_status)
    disp(append("Warning: end point not reached for ", strjoin(date_label(~channel_status), ", "), ...
        " (continued with the centerline of the last attempt; check these centerlines, or try other search parameters). "))
end


%% all dates together
%  on the most recent DEM (still in memory), one color per date (dark
%  purple: first date, yellow: last date)

% geometries in image coordinates of the most recent DEM
x_cent_r = cell(no_DEMs, 1);  y_cent_r = cell(no_DEMs, 1);
x_prof_r = cell(no_DEMs, 1);  y_prof_r = cell(no_DEMs, 1);
edge_coord_r = cell(no_DEMs, 1);
for t = 1:no_DEMs
    [x_cent_r{t}, y_cent_r{t}] = worldToIntrinsic(R, xw_cent{t}, yw_cent{t});
    [x_prof_r{t}, y_prof_r{t}] = worldToIntrinsic(R, xw_prof{t}, yw_prof{t});
    edge_coord_r{t} = NaN(size(edge_coord_w{t}));
    [edge_coord_r{t}(:,1), edge_coord_r{t}(:,2)] = worldToIntrinsic(R, edge_coord_w{t}(:,1), edge_coord_w{t}(:,2));
    [edge_coord_r{t}(:,3), edge_coord_r{t}(:,4)] = worldToIntrinsic(R, edge_coord_w{t}(:,3), edge_coord_w{t}(:,4));
end
date_cols = sty.cmap_prof_inferno(no_DEMs);

f_all_map = figure;
plot_overview_map(axes(f_all_map), DEM, R, x_cent_r, y_cent_r, edge_coord_r, x_prof_r, y_prof_r, date_label, ...
    'clims', map_clims, 'hillshade', map_hillshade, 'colors', date_cols, ...
    'bridge_gaps', plot_edge_gaps == 1, 'title', "Mapped channel, all dates", ...
    'subtitle', append(string(no_DEMs), " dates, on the DEM of ", date_label(end), ", edge detection: ", edge_method_txt))

f_all = new_figure([18 16], save_figs);
plot_along_all_figure(f_all, T_along, cellfun(@(P) mean(P, 1, 'omitnan')', profiles, 'UniformOutput', false), ...
    date_label, 'colors', date_cols, 'title', "All dates", 'subtitle', "Mean profile elevation, depth and width along the channel")


%% write to files

if save_figs || save_shps || save_tables
    disp("Writing figures, shapefiles and/or tables... ")

    if save_figs
        print(f_all_map, append(fig_dir, file_prefix, 'all_dates_mapped_channel'), figs_filetype, figs_resolution)
        print(f_all, append(fig_dir, file_prefix, 'all_dates_along_channel'), figs_filetype, figs_resolution)
        close(f_all)
    end

    % write mapped geometries to shapefile (in image coordinates of the
    % most recent DEM, i.e. its spatial referencing R)
    if save_shps

        % all centerlines in a single file
        % attributes: date, start point attributes (if any), whether the
        % end point was reached, and the search parameters actually used
        % (those of the successful attempt, incl. retries, see find_centerline())
        % (for loaded centerlines: the attributes of the input shapefile)
        attr = table(date_label(:), string(DEM_files), 'VariableNames', {'date', 'dem_file'});
        attr = [attr, attr_startend];
        if start_end_method ~= 4
            info = struct2table(vertcat(cent_info{:}), 'AsArray', true);
            info = renamevars(info, info.Properties.VariableNames, ...
                ["step_m", "angle_deg", "max_grad", "win_cent", "reversed", "perturbed", "attempt", "end_found"]);
            info.win_dem = repmat(window_DEM, no_DEMs, 1);
            attr = [attr, info];
        end
        fn = append(shp_dir, file_prefix, 'all_centerlines');
        lines_to_shp(x_cent_r, y_cent_r, R, fn, 'channel_label', num2cell(date_label), 'attributes', attr);

        % centerline and outlines in one file per date
        for t = 1:no_DEMs
            outlines_x = {x_cent_r{t}; edge_coord_r{t}(:,1); edge_coord_r{t}(:,3)};
            outlines_y = {y_cent_r{t}; edge_coord_r{t}(:,2); edge_coord_r{t}(:,4)};
            fline = {"centerline", "left_edge", "right_edge"};
            fn = append(shp_dir, file_prefix, date_label(t), "_outlines");
            lines_to_shp(outlines_x, outlines_y, R, fn, 'line_type', fline);
        end

        % all profile transects in one file per date
        for t = 1:no_DEMs
            fn = append(shp_dir, file_prefix, date_label(t), "_profiles");
            lines_to_shp(x_prof_r{t}, y_prof_r{t}, R, fn, 'prof_no', 1:size(x_prof_r{t}, 2));
        end
    end

    % write profiles and channel metrics to tables (.csv)
    if save_tables
        T_sum = table();
        for t = 1:no_DEMs
            % channel metrics per profile (see channel_metrics())
            writetable(T_along{t}, append(tab_dir, file_prefix, date_label(t), '_metrics.csv'))

            % cross sectional profile elevations: one column per profile,
            % first column: distance from profile center [m] (positive =
            % left side, looking downstream)
            prof_dist = ((1:size(profiles{t}, 1))' - (size(profiles{t}, 1)+1)/2)*res_t(t);
            T_prof = array2table([prof_dist, profiles{t}], 'VariableNames', ...
                ["dist_center_m", append("prof_", string(1:size(profiles{t}, 2)))]);
            writetable(T_prof, append(tab_dir, file_prefix, date_label(t), '_profiles.csv'))

            % summary statistics per date
            T_sum = [T_sum; channel_summary(T_along{t}, date_label(t), channel_status(t), channel_length{t})];
        end
        T_sum = renamevars(T_sum, 'channel', 'date');
        writetable(T_sum, append(tab_dir, file_prefix, 'summary_per_date.csv'))
    end

    disp(append("Done writing files. Check '", append(results_dir, proj_subdir), "' for output. "))
end

disp("Done! ")
