clear all
close all
clc

%% Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
%  Algorithm that returns ice shelf basal channels' centerline, outlines
%  and cross sectional profiles; based on surface expressions in a DEM. 
%  This script maps one or more channels given a single DEM and user
%  identified channel start and end points, or existing channel centerlines
%  (e.g. selected from a previous run). (See "map_channel_timeseries.m" for
%  mapping a single channel over different DEMs.)
% 
%  Ensemble mode: if sets of values are given for any of the centerline
%  search parameters in config.m, centerlines are mapped for every
%  combination of values (full-factorial ensemble). Only centerlines that
%  connect start and end point are kept, and written to a single shapefile
%  incl. their search parameters. The script then stops: select a
%  centerline per channel (e.g. by visual inspection) and run again with
%  start_end_method = 4 to extract profiles and edges.
% 
%  Input: 
%   - DEM (GeoTIFF) 
%   - channel start and end points (through GUI or read from shapefile), or
%     channel centerlines (shapefile)
% 
%  Output: 
%   - channel centerlines and outlines (figures and georeferenced shapefiles)
%     (ensemble mode: centerline ensemble only)
%   - cross sectional elevation/depth profiles (figures and tables)
%   - channel depth and width along centerlines (figures and tables)
%   - summary of channel metrics (table)
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
disp(append("Running Sub-Ice with configuration file ", config_fn, ". "))
% (evaluated in the current folder: relative paths in the configuration
% file are relative to the Sub-Ice folder, wherever the file is stored)
eval(fileread(config_fn))

% ensemble mode if sets of values are given for any of the search parameters
ensemble_mode = any(cellfun(@numel, {window_DEM, search_step, search_angle, max_gradient, window_cent}) > 1); 
if ensemble_mode && start_end_method == 4
    error("Ensemble mode (sets of centerline search parameters in config.m) requires a centerline search: set start_end_method to 1, 2 or 3. ")
end


%% read DEM from GeoTIFF
%  and smooth if required

if path_to_DEM(end-3:end) ~= ".tif"
    error("path_to_DEM should point to a single DEM .tif file (check in config.m). ")
end

% read from geotiff
[DEM, R] = readgeoraster(path_to_DEM);
res = R.CellExtentInWorldX;     % resolution of DEM [m]

% replace no data with NaN
DEM(DEM==DEM_nodata) = NaN;
disp("Loaded DEM. ")

% smooth DEM for centerline mapping, if requested
% (smoothed DEM is cached to disk and loaded from there in subsequent runs;
% in ensemble mode: for every value of window_DEM, see below)
if ~ensemble_mode && start_end_method ~= 4
    DEM_smooth = smooth_dem(DEM, res, window_DEM, path_to_DEM);
end


%% specify channel centerline start/end points
%  have user click on channel start and end points through GUI, read start
%  and end points from shapefile or enter img coordinates manually; or load
%  existing channel centerlines from shapefile

% figure style (colors, font sizes; see fig_style())
sty = fig_style(); 

% map of the DEM, to show start/end points (see plot_map_base())
figure(1)
plot_map_base(gca, DEM, R, 'clims', map_clims, 'hillshade', map_hillshade); 
text_offs = 3;  % [pix]

if start_end_method == 1 
    % click on start/end points 
    c = 0; 
    while true
        c = c+1; 
        title(append("Left click on channel start and end points. Right click when done! "))
        disp(append("Left click on channel start and end points. Right click when done! "))
        
        % start point
        [P_start(c,1), P_start(c,2), button] = ginput(1);
        if button == 3                  % right click
            c = c-1; 
            P_start = P_start(1:c, :);  % keep valid entries only
            break
        end
        scatter(P_start(c,1), P_start(c,2), 'm', 'filled')
        text(P_start(c,1) + text_offs, P_start(c,2) + text_offs, 'start', 'Color', 'm')        
        
        % end point
        [P_end(c,1), P_end(c,2), button] = ginput(1); 
        if button == 3
            c = c-1; 
            P_start = P_start(1:c, :);
            P_end = P_end(1:c, :);
            break
        end
        text(P_end(c,1) + text_offs, P_end(c,2) + text_offs, 'end', 'Color', 'm')
        scatter(P_end(c,1), P_end(c,2), 'm', 'filled')
    end
    
    no_channels = c; 
    disp(append("Entered ", string(no_channels), " channels' start and end points. "))
    channel_label = strings(no_channels, 1);
    for c = 1:no_channels
        channel_label(c) = string(input(append("Provide a label for channel ", string(c), ': '), 's'));
    end
    
    
elseif start_end_method == 2
    % read start/end points from shapefile
    % only works if shapefile has same map projection as DEM!
    S = shaperead(path_to_start_end_shp); 
    % filter on "shelf" field
    if ~isempty(shelf_filter)
        S = S(strcmp({S.Shelf}, shelf_filter));
    end
    [x_startend, y_startend] =  worldToIntrinsic(R, vertcat(S.X), vertcat(S.Y));
    
    % assumption: points are ordered, start channel 1, end channel 1, ... 
    % > shapefile should have even amount of points
    if rem(length(x_startend), 2) == 0
        idx = 1:length(x_startend); 
        even_idx = rem(idx, 2) == 0; 
        P_start = [x_startend(~even_idx), y_startend(~even_idx)]; 
        P_end = [x_startend(even_idx), y_startend(even_idx)]; 
    else
        error('Found odd number of start and end points in shapefile. Please provide an even number of points, ordered channel start, end respectively. '); 
    end 

    scatter(P_start(:,1), P_start(:,2), 'm', 'filled')
    scatter(P_end(:,1), P_end(:,2), 'm', 'filled')
    text(P_start(:,1) + text_offs, P_start(:,2) + text_offs, 'start', 'Color', 'm')
    text(P_end(:,1) + text_offs, P_end(:,2) + text_offs, 'end', 'Color', 'm')

    no_channels = size(P_start, 1); 
    disp(append("Found ", string(no_channels), " channels' start and end points in shapefile. "))
    
    % channel label - feel free to adjust: 
    channel_label = append("channel-", string(1:no_channels)); 
    % (should be string)

    % attributes of start points (e.g. qualitative flags added by the
    % user), carried over to the centerline shapefile
    % (start and end point of a channel are assumed to share attributes)
    attr_startend = struct2table(S(~even_idx), 'AsArray', true);
    attr_startend = removevars(attr_startend, intersect(attr_startend.Properties.VariableNames, {'Geometry', 'BoundingBox', 'X', 'Y'}));
    
    
elseif start_end_method == 3
    % enter start/end points manually below
    P_start = [353, 128];    % x, y in image coord [pix]
    P_end = [260, 84];       % x, y in image coord [pix]
    scatter(P_start(:,1), P_start(:,2), 'm', 'filled')
    scatter(P_end(:,1), P_end(:,2), 'm', 'filled')
    text(P_start(:,1) + text_offs, P_start(:,2) + text_offs, 'start', 'Color', 'm')
    text(P_end(:,1) + text_offs, P_end(:,2) + text_offs, 'end', 'Color', 'm')

    no_channels = size(P_start, 1); 
    channel_label = append("channel-", string(1:no_channels)); 


elseif start_end_method == 4
    % load existing centerlines from line shapefile (no centerline search)
    % only works if shapefile has same map projection as DEM!
    S = shaperead(path_to_centerline_shp); 
    no_channels = length(S); 

    x_load = cell(no_channels, 1); 
    y_load = cell(no_channels, 1); 
    for c = 1:no_channels
        % strip NaN part separators and convert to img coordinates
        if sum(isnan(S(c).X(1:end-1))) > 0
            disp(append("Warning: centerline ", string(c), " in shapefile consists of multiple parts. Parts are joined. "))
        end
        valid = ~isnan(S(c).X) & ~isnan(S(c).Y); 
        [x_load{c}, y_load{c}] = worldToIntrinsic(R, S(c).X(valid)', S(c).Y(valid)'); 
    end
    P_start = cell2mat(cellfun(@(x, y) [x(1), y(1)], x_load, y_load, 'UniformOutput', false)); 
    P_end = cell2mat(cellfun(@(x, y) [x(end), y(end)], x_load, y_load, 'UniformOutput', false)); 

    % channel label from attribute table, if available (first attribute
    % starting with "channel"; field names are truncated on write)
    attr_startend = struct2table(S, 'AsArray', true);
    attr_startend = removevars(attr_startend, intersect(attr_startend.Properties.VariableNames, {'Geometry', 'BoundingBox', 'X', 'Y'}));
    label_field = attr_startend.Properties.VariableNames(startsWith(attr_startend.Properties.VariableNames, "channel", 'IgnoreCase', true)); 
    if isempty(label_field)
        channel_label = append("channel-", string(1:no_channels)); 
    else
        channel_label = string(attr_startend.(label_field{1}))'; 
        attr_startend = removevars(attr_startend, label_field{1}); % (written separately)
    end
    % (other attributes, e.g. search parameters, are carried over to the
    % output shapefile)

    scatter(P_start(:,1), P_start(:,2), 'm', 'filled')
    scatter(P_end(:,1), P_end(:,2), 'm', 'filled')
    text(P_start(:,1) + text_offs, P_start(:,2) + text_offs, 'start', 'Color', 'm')
    text(P_end(:,1) + text_offs, P_end(:,2) + text_offs, 'end', 'Color', 'm')

    disp(append("Loaded ", string(no_channels), " channel centerlines from shapefile. "))


else
    error('Specify a valid method to enter channel start/end points (under "user specifiable variables" in config.m). ')
end

% no start point attributes, unless read from shapefile
if ~exist('attr_startend', 'var')
    attr_startend = table('Size', [no_channels 0], 'VariableTypes', {});
end

% exclude channels, if requested
if ~isempty(excl_channel_labels)
    keep = ~ismember(channel_label, excl_channel_labels); 
    disp(append("Excluding ", string(sum(~keep)), " channel(s): ", strjoin(channel_label(~keep), ", "), ". "))
    P_start = P_start(keep, :); 
    P_end = P_end(keep, :); 
    channel_label = channel_label(keep); 
    attr_startend = attr_startend(keep, :); 
    if start_end_method == 4
        x_load = x_load(keep); 
        y_load = y_load(keep); 
    end
    no_channels = sum(keep); 
end


%% ensemble mode: map centerline ensemble
%  one centerline per channel for every combination of search parameter
%  values; keep those that connect start and end point; write to file; stop

if ensemble_mode

    % all combinations of search parameter values (full-factorial ensemble)
    [g_step, g_angle, g_grad, g_wcent] = ndgrid(search_step, search_angle, max_gradient, window_cent); 
    no_comb = numel(g_step); 
    no_members = no_comb*numel(window_DEM); 
    disp(append("Ensemble mode: mapping ", string(no_members), " centerline realizations for each of ", string(no_channels), " channels. "))
    if ~save_shps
        disp("Warning: save_shps = 0, the centerline ensemble will not be written to file. ")
    end

    % store realizations that connect start and end point
    ens_x = {}; 
    ens_y = {}; 
    ens_channel = []; 
    ens_info = struct([]); 
    ens_set = []; 
    no_success = zeros(no_channels, 1); 

    member = 0; 
    for w = 1:numel(window_DEM)
        DEM_smooth = smooth_dem(DEM, res, window_DEM(w), path_to_DEM);
        for k = 1:no_comb
            member = member + 1; 
            disp(append("Ensemble member ", string(member), "/", string(no_members), ...
                ": window_DEM = ", string(window_DEM(w)), ", search_step = ", string(g_step(k)), ", search_angle = ", string(g_angle(k)), ...
                ", max_gradient = ", string(g_grad(k)), ", window_cent = ", string(g_wcent(k)), ". "))
            for c = 1:no_channels
                [x_c, y_c, ~, info_c] = find_centerline(P_start(c,:), P_end(c,:), DEM_smooth, R, ... 
                                                    'search_step',      g_step(k), ... 
                                                    'search_angle',     g_angle(k), ... 
                                                    'max_gradient',     g_grad(k), ... 
                                                    'window',           g_wcent(k), ...
                                                    'max_recursions',   max_recursions); 
                if info_c.end_found
                    ens_x{end+1, 1} = x_c; 
                    ens_y{end+1, 1} = y_c; 
                    ens_channel(end+1, 1) = c; 
                    info_c.win_dem = window_DEM(w); 
                    ens_info = [ens_info; info_c]; 
                    % configured (not perturbed) parameter values and member no.
                    ens_set(end+1, :) = [member, window_DEM(w), g_step(k), g_angle(k), g_grad(k), g_wcent(k)]; 
                    no_success(c) = no_success(c) + 1; 
                end
            end
        end
    end

    % success rate per channel
    disp("Ensemble realizations connecting start and end point, per channel: ")
    T_success = table(channel_label(:), no_success, repmat(no_members, no_channels, 1), round(100*no_success/no_members, 1), ...
        'VariableNames', {'channel', 'no_realizations', 'no_members', 'success_rate_pct'}); 
    disp(T_success)
    if save_tables
        writetable(T_success, append(tab_dir, file_prefix, 'centerline_ensemble_success.csv'))
    end

    % overview figure: all realizations, one color per channel, on the DEM
    % (map extent: realizations and start/end points, plus a margin)
    f = figure(1); 
    clf(f)
    ax_map = axes(f); 
    all_x = [vertcat(ens_x{:}); P_start(:,1); P_end(:,1)]; 
    all_y = [vertcat(ens_y{:}); P_start(:,2); P_end(:,2)]; 
    margin = 0.1*max(max(all_x) - min(all_x), max(all_y) - min(all_y)); 
    plot_map_base(ax_map, DEM, R, 'clims', map_clims, 'hillshade', map_hillshade, ...
        'extent', [min(all_x)-margin, max(all_x)+margin, min(all_y)-margin, max(all_y)+margin]); 
    ch_cols = sty.cmap_channels(no_channels); 
    for e = 1:numel(ens_x)
        plot(ax_map, ens_x{e}, ens_y{e}, 'Color', [ch_cols(ens_channel(e),:) 0.6], 'LineWidth', sty.lw_thin, 'HandleVisibility', 'off')
    end
    scatter(ax_map, [P_start(:,1); P_end(:,1)], [P_start(:,2); P_end(:,2)], 12, sty.ink, 'filled', 'HandleVisibility', 'off')
    for c = 1:no_channels
        text(ax_map, P_start(c,1), P_start(c,2), append(" ", channel_label(c), " "), 'Color', sty.label, ...
            'FontSize', sty.font_small, 'BackgroundColor', 'w', 'Margin', 0.5, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
    panel_title(ax_map, "Centerline ensemble", append(string(numel(ens_x)), " realizations connecting start and end point (", ...
        string(no_members), " parameter combinations x ", string(no_channels), " channels)"))

    if save_figs
        fn = append(fig_dir, file_prefix, 'centerline_ensemble'); 
        print(fn, figs_filetype, figs_resolution)
    end

    % write ensemble to shapefile: all realizations of all channels in one
    % file, ordered by channel; attributes: start point attributes (if any),
    % search parameters actually used (incl. retries, see find_centerline()),
    % ensemble member no. and the configured parameter values of the member
    if save_shps && ~isempty(ens_x)
        [~, order] = sort(ens_channel); 
        info = struct2table(ens_info(order), 'AsArray', true);
        info = renamevars(info, ["search_step", "search_angle", "max_gradient", "window"], ...
            ["step_m", "angle_deg", "max_grad", "win_cent"]);
        set_info = array2table(ens_set(order, :), 'VariableNames', ...
            ["member", "set_wdem", "set_step", "set_angle", "set_grad", "set_wcent"]); 
        attr = [attr_startend(ens_channel(order), :), info, set_info];
        fn = append(shp_dir, file_prefix, 'centerline_ensemble'); 
        lines_to_shp(ens_x(order), ens_y(order), R, fn, 'channel_label', num2cell(channel_label(ens_channel(order))), 'attributes', attr); 
        disp(append("Wrote centerline ensemble to ", fn, ".shp. "))
    elseif save_shps
        disp("No realization connected start and end point, no shapefile written. ")
    end

    disp("Done! Select a centerline per channel from the ensemble and run again with start_end_method = 4 to extract profiles and channel edges. ")
    return
end


%% map channel geometries and extract profiles
%  and visualize on overview figure

x_cent = cell(no_channels, 1);  
y_cent = cell(no_channels, 1); 
channel_length = cell(no_channels, 1); 
profiles = cell(no_channels, 1); 
x_prof = cell(no_channels, 1); 
y_prof = cell(no_channels, 1); 
edge_idx = cell(no_channels, 1); 
edge_coord = cell(no_channels, 1); 
edge_elev = cell(no_channels, 1); 
valid_edges = cell(no_channels, 1); 
fchannel = cell(no_channels, 1);
cent_info = cell(no_channels, 1); 
trough_elev = cell(no_channels, 1); 
trough_depth = cell(no_channels, 1); 
valid_trough = cell(no_channels, 1); 

% to store whether we "successfully" found the channel centerline: 
channel_status = zeros(no_channels, 1); 

% loop over channels
for c = 1:no_channels
    disp(append("Start mapping channel geometry of ", channel_label(c), ". ")) 
    
    if start_end_method == 4
        % centerline loaded from shapefile
        x_cent{c} = x_load{c}; 
        y_cent{c} = y_load{c}; 
        cent_length = sqrt(diff(x_cent{c}).^2 + diff(y_cent{c}).^2); 
        channel_status(c) = 1; 

    else
        % find centerline
        % attention: possibly using smoothed DEM to find centerlines! 
        % (cent_info: search parameters of the (last) attempt, incl. retries)
        [x_cent{c}, y_cent{c}, cent_length, cent_info{c}] = find_centerline(P_start(c,:), P_end(c,:), DEM_smooth, R, ... 
                                                    'search_step',      search_step, ... 
                                                    'search_angle',     search_angle, ... 
                                                    'max_gradient',     max_gradient, ... 
                                                    'window',           window_cent, ...
                                                    'max_recursions',   max_recursions); 

        % did we reach the channel end point?
        % (if not, continue with the centerline of the last attempt)
        channel_status(c) = cent_info{c}.end_found;
        if ~channel_status(c)
            disp(append("Warning: end point of ", channel_label(c), " not reached. Continuing with the centerline of the last attempt (", string(round(sum(cent_length)*res/1000, 1)), " km). "))
        end
    end
    channel_length{c} = sum(cent_length);       % in [pix]
    channel_length{c} = channel_length{c}*res;  % in [m]

    % find cross sectional profiles
    [profiles{c}, x_prof{c}, y_prof{c}] = find_profiles(x_cent{c}, y_cent{c}, DEM, R, ...
                                                'prof_length',      prof_length, ...
                                                'prof_interval',    prof_interval);
    no_profiles = size(profiles, 2); 

    % find channel edges/outlines
    [edge_idx{c}, edge_coord{c}, edge_elev{c}, valid_edges{c}] = find_edges(profiles{c}, x_prof{c}, y_prof{c}, res, ...
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
                                                'edge_subst_window', edge_subst_window);

    % find channel trough elevation and depth
    [trough_elev{c}, trough_depth{c}, valid_trough{c}] = find_trough(profiles{c}, edge_elev{c}, ...
                                                'z_thr',            z_thr_trough_elev, ...
                                                'subst_window',     trough_subst_window); 

    % label for shapefile
    fchannel{c} = channel_label(c); 
    
end

disp(append("Finished mapping. End point reached for ", string(sum(channel_status)), "/", string(length(channel_status)), " channels. "))


%% overview figure
%  channel centerlines and the area between the channel edges, on the DEM
%  (map extent: all profiles, plus a margin)

f = figure(1); 
clf(f)
ax_map = axes(f); 
all_x = [x_prof{:}]; 
all_y = [y_prof{:}]; 
all_x = all_x([1 end], :); 
all_y = all_y([1 end], :); 
margin = 0.05*max(max(all_x(:)) - min(all_x(:)), max(all_y(:)) - min(all_y(:))); 
plot_map_base(ax_map, DEM, R, 'clims', map_clims, 'hillshade', map_hillshade, ...
    'extent', [min(all_x(:))-margin, max(all_x(:))+margin, min(all_y(:))-margin, max(all_y(:))+margin]); 

for c = 1:no_channels
    % profile transects (thin line from first to last profile point)
    if plot_prof_transects == 1
        plot(ax_map, x_prof{c}([1 end], :), y_prof{c}([1 end], :), 'Color', [sty.ink 0.25], 'LineWidth', 0.5, 'HandleVisibility', 'off')
    end
    % channel area between the edges (with outlines along the edges)
    shade_channel(ax_map, edge_coord{c}, sty.shade, sty.shade_alpha, sty.shade_outline*sty.shade, plot_edge_gaps == 1)
end
for c = 1:no_channels
    % centerline
    plot(ax_map, x_cent{c}, y_cent{c}, 'Color', sty.cent, 'LineWidth', sty.lw_data, 'HandleVisibility', 'off')
end
for c = 1:no_channels
    % channel label at the upstream end
    text(ax_map, x_cent{c}(1), y_cent{c}(1), append(" ", channel_label(c), " "), 'Color', sty.label, ...
        'FontSize', sty.font_small, 'BackgroundColor', 'w', 'Margin', 0.5, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end

% legend (proxy entries)
h_leg = [plot(ax_map, NaN, NaN, 'Color', sty.cent, 'LineWidth', sty.lw_data), ...
         patch(ax_map, NaN, NaN, sty.shade, 'FaceAlpha', sty.shade_alpha, 'EdgeColor', sty.shade_outline*sty.shade)]; 
leg_txt = ["Centerline", "Channel (between edges)"]; 
if plot_prof_transects == 1
    h_leg(end+1) = plot(ax_map, NaN, NaN, 'Color', [sty.ink 0.25], 'LineWidth', 0.5); 
    leg_txt(end+1) = "Cross-sectional profiles"; 
end
legend(ax_map, h_leg, leg_txt, 'Location', 'northeast', 'FontSize', sty.font_small, 'TextColor', sty.ink, 'Box', 'on', 'EdgeColor', sty.grid)
edge_method_txt = string(edge_method); 
if edge_method ~= "SlopeThreshold"
    edge_method_txt = append(edge_method_txt, "-", knee_method);   % (knee point method also used by NearPeaks)
end
panel_title(ax_map, "Mapped channels", append(string(no_channels), " channels, edge detection: ", edge_method_txt))


%% write to files

if save_figs || save_shps || save_tables
    disp("Writing figures, shapefiles and/or tables... ")

    % print overview figure to file
    if save_figs
        %f.WindowState = 'maximized'; % make figure fullscreen before saving
        fn = append(fig_dir, file_prefix, 'mapped_channels'); 
        print(fn, figs_filetype, figs_resolution)
        %f.WindowState = 'normal'; 
    end
    
    % write mapped geometries to shapefile
    if save_shps
    
        % all centerlines in a single file
        % attributes: start point attributes (if any), whether the end
        % point was reached, and the search parameters actually used
        % (those of the successful attempt, incl. retries, see find_centerline())
        % (for loaded centerlines: the attributes of the input shapefile)
        if start_end_method == 4
            attr = attr_startend; 
        else
            info = struct2table(vertcat(cent_info{:}), 'AsArray', true);
            info = renamevars(info, info.Properties.VariableNames, ...
                ["step_m", "angle_deg", "max_grad", "win_cent", "reversed", "perturbed", "attempt", "end_found"]);
            info.win_dem = repmat(window_DEM, no_channels, 1);
            attr = [attr_startend, info];
        end
        fn = append(shp_dir, file_prefix, 'all_centerlines'); 
        lines_to_shp(x_cent, y_cent, R, fn, 'channel_label', fchannel, 'attributes', attr); 
        
        % centerlines and outlines in one file per channel
        for c = 1:no_channels
            outlines_x = cell(3, 1); 
            outlines_y = cell(3, 1); 
            outlines_x{1} = x_cent{c};            % centerline
            outlines_y{1} = y_cent{c}; 
            outlines_x{2} = edge_coord{c}(:,1);   % left edge
            outlines_y{2} = edge_coord{c}(:,2); 
            outlines_x{3} = edge_coord{c}(:,3);   % right edge
            outlines_y{3} = edge_coord{c}(:,4); 
            fline = {"centerline", "left_edge", "right_edge"}; 
            fn = append(shp_dir, file_prefix, channel_label(c), "_outlines"); 
            lines_to_shp(outlines_x, outlines_y, R, fn, 'line_type', fline);
        end
        
        % all profile transects in one file per channel
        for c = 1:no_channels
            no_profiles = size(x_prof{c}, 2); 
            fprof = 1:no_profiles;  
            fn = append(shp_dir, file_prefix, channel_label(c), "_profiles"); 
            lines_to_shp(x_prof{c}, y_prof{c}, R, fn, 'prof_no', fprof);
        end 
    end

    % write profiles and channel metrics to tables (.csv)
    if save_tables
        T_sum = table(); 
        for c = 1:no_channels
            % channel metrics per profile (see channel_metrics())
            T = channel_metrics(x_cent{c}, y_cent{c}, x_prof{c}, y_prof{c}, ...
                edge_idx{c}, edge_coord{c}, edge_elev{c}, valid_edges{c}, ...
                trough_elev{c}, valid_trough{c}, R, prof_interval); 
            writetable(T, append(tab_dir, file_prefix, channel_label(c), '_metrics.csv'))

            % cross sectional profile elevations: one column per profile,
            % first column: distance from profile center [m] (positive =
            % left side, looking downstream)
            prof_dist = ((1:size(profiles{c}, 1))' - (size(profiles{c}, 1)+1)/2)*res; 
            T_prof = array2table([prof_dist, profiles{c}], 'VariableNames', ...
                ["dist_center_m", append("prof_", string(1:size(profiles{c}, 2)))]); 
            writetable(T_prof, append(tab_dir, file_prefix, channel_label(c), '_profiles.csv'))

            % summary statistics per channel
            T_sum = [T_sum; table(channel_label(c), channel_status(c), channel_length{c}/1000, height(T), ...
                mean(T.width_m, 'omitnan'), median(T.width_m, 'omitnan'), min(T.width_m), max(T.width_m), std(T.width_m, 'omitnan'), ...
                mean(T.hw_left_m, 'omitnan'), mean(T.hw_right_m, 'omitnan'), ...
                mean(T.depth_m, 'omitnan'), median(T.depth_m, 'omitnan'), max(T.depth_m), ...
                mean(T.elev_ledge_m, 'omitnan'), mean(T.elev_redge_m, 'omitnan'), mean(T.elev_trough_m, 'omitnan'), ...
                100*mean(~T.valid_ledge), 100*mean(~T.valid_redge), 100*mean(~T.valid_trough), ...
                'VariableNames', {'channel', 'end_found', 'length_km', 'no_profiles', ...
                'mean_width_m', 'median_width_m', 'min_width_m', 'max_width_m', 'std_width_m', ...
                'mean_hw_left_m', 'mean_hw_right_m', ...
                'mean_depth_m', 'median_depth_m', 'max_depth_m', ...
                'mean_elev_ledge_m', 'mean_elev_redge_m', 'mean_elev_trough_m', ...
                'outliers_ledge_pct', 'outliers_redge_pct', 'outliers_trough_pct'})]; 
        end
        writetable(T_sum, append(tab_dir, file_prefix, 'channel_summary.csv'))
    end

    disp(append("Done writing files. Check '", append(results_dir, proj_subdir), "' for output. "))
end


%% extended figures
%  cross sectional profiles, metrics along channel length etc.

if ext_figs
disp("Creating and possibly saving extended figures. Sit tight. ")
    T_along = cell(no_channels, 1);     % channel metrics per profile, per channel
    
    for c = 1:no_channels
        
        % cross sectional profiles
        no_profiles = size(profiles{c}, 2); 
        
        % distance of profile samples from profile center [m] (positive = left)
        prof_dist_vector = ((1:size(profiles{c}, 1)) - (size(profiles{c}, 1)+1)/2)*res; 

        % full cross sectional profiles, absolute elevation
        figure('Units', 'centimeters', 'Position', [2 2 18 10])
        plot_profile_figure(gca, prof_dist_vector, profiles{c}, edge_idx{c}, edge_elev{c}, ...
            'title', channel_label(c), ...
            'subtitle', append("Cross-sectional profiles (", string(no_profiles), "), looking upstream"), ...
            'min_width', min_width, 'max_width', max_width)

        if save_figs
            fn = append(fig_dir, file_prefix, channel_label(c), '_full_profiles_elev');
            print(fn, figs_filetype, figs_resolution)
        end

        % cross sectional profiles between the channel edges, relative to
        % the left channel edge
        figure('Units', 'centimeters', 'Position', [2 2 18 10])
        plot_profile_figure(gca, prof_dist_vector, profiles{c}, edge_idx{c}, edge_elev{c}, ...
            'relative', true, 'title', channel_label(c), ...
            'subtitle', append("Cross-sectional profiles between the channel edges, looking upstream; mean depth below left edge: ", ...
                string(round(mean(edge_elev{c}(:,1) - trough_elev{c}, 'omitnan'), 1)), " m"), ...
            'min_width', min_width, 'max_width', max_width)

        if save_figs
            fn = append(fig_dir, file_prefix, channel_label(c), '_lim_profiles_depth'); 
            print(fn, figs_filetype, figs_resolution)
        end
        
        
        % channel geometry along the channel (one figure per channel)
        % (metrics per profile, see channel_metrics())
        T_along{c} = channel_metrics(x_cent{c}, y_cent{c}, x_prof{c}, y_prof{c}, ...
            edge_idx{c}, edge_coord{c}, edge_elev{c}, valid_edges{c}, ...
            trough_elev{c}, valid_trough{c}, R, prof_interval); 
        f_along = figure('Units', 'centimeters', 'Position', [2 2 18 13]); 
        plot_along_channel_figure(f_along, T_along{c}, channel_label(c))

        if save_figs
            fn = append(fig_dir, file_prefix, channel_label(c), '_along_channel');
            print(f_along, fn, figs_filetype, figs_resolution)
        end
    end
    
    % all channels together: mean profile elevation, depth and width along
    % the channel (one color per channel)
    f_all = figure('Units', 'centimeters', 'Position', [2 2 18 16]); 
    tl = tiledlayout(f_all, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact'); 
    ch_cols = sty.cmap_channels(no_channels); 
    ax_all = gobjects(3, 1); 
    y_lab = {'Mean profile elevation [m]', 'Depth below left edge [m]', 'Width [km]'}; 
    for k = 1:3
        ax_all(k) = nexttile(tl); 
        hold(ax_all(k), 'on')
        for c = 1:no_channels
            d = T_along{c}.dist_along_m/1000; 
            switch k
                case 1, y = mean(profiles{c}, 1, 'omitnan')'; 
                case 2, y = T_along{c}.depth_m; 
                case 3, y = T_along{c}.width_m/1000; 
            end
            plot(ax_all(k), d, y, 'Color', ch_cols(c,:), 'LineWidth', sty.lw_data)
        end
        style_axis(ax_all(k))
        ylabel(ax_all(k), y_lab{k}, 'Color', sty.ink)
        if k == 2
            ax_all(k).YDir = 'reverse';     % (depth positive down: deeper channel = lower line)
        end
        if k < 3
            xticklabels(ax_all(k), {})
        end
    end
    xlabel(ax_all(3), 'Distance along channel [km]', 'Color', sty.ink)
    linkaxes(ax_all, 'x')
    xlim(ax_all(1), [0, max(cellfun(@(T) max([T.dist_along_m; 1]), T_along))/1000])
    panel_title(ax_all(1), "All channels", "Mean profile elevation, depth and width along the channel")
    legend(ax_all(1), channel_label, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
        'Box', 'off', 'FontSize', sty.font_small, 'TextColor', sty.ink)
    
    if save_figs
        fn = append(fig_dir, file_prefix, 'all_channels_along_channel'); 
        print(f_all, fn, figs_filetype, figs_resolution)
    end
    
end
    
disp("Done! ")
