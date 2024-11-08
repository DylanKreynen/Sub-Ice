clear all
close all
clc

%% Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
%  Algorithm that returns ice shelf basal channels' centerline, outlines
%  and cross sectional profiles; based on surface expressions in a DEM. 
%  This script maps one or more channels given a single DEM and user
%  identified channel start and end points. (See "map_channel_timeseries.m"
%  for mapping a single channel over different DEMs.)
% 
%  Input: 
%   - DEM (GeoTIFF) 
%   - channel start and end points (through GUI or read from shapefile)
% 
%  Output: 
%   - channel centerlines and outlines (figures and georeferenced shapefiles)
%   - cross sectional elevation/depth profiles (figures)
%   - channel depth and width along centerlines (figures)
% 
%  (c) Dylan Kreynen
%  University of Oslo
%  2024
% 
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% run configuration file - set user specifiable parameters
%  contains filepaths to input/output data, sets method to select channel
%  start/end points, configures centerline search parameters, output 
%  behaviour etc. > update configuration file as required! 

run 'config.m'


%% read DEM from GeoTIFF
%  and some basic manipulation

[DEM, R] = readgeoraster(path_to_DEM);
res = R.CellExtentInWorldX;     % resolution of DEM [m]

% remove no data values, aid visualization
% (update as required)
DEM(DEM<-20) = NaN; % no data: -9999
DEM(DEM>100) = 50; 


%% specify channel centerline start/end points
%  have user click on channel start and end points through GUI, read start
%  and end points from shapefile or enter img coordinates manually

figure(1)
imagesc(DEM)
hold on
axis image
colormap gray
xlabel('x coord [pix]')
ylabel('y coord [pix]')
text_offs = 3;  % [pix]

if start_end_method == 1 
    % click on start/end points 
    c = 0; 
    while 1 == 1
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
    channel_label = strings(no_channels);
    for c = 1:no_channels
        channel_label(c) = string(input(append("Provide a label for channel ", string(c), ': '), 's'));
    end
    
    
elseif start_end_method == 2
    % read start/end points from shapefile
    % only works if shapefile has same map projection as DEM!
    S = shaperead(path_to_start_end_shp); 
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
    
    
elseif start_end_method == 3
    % enter start/end points manually below
    P_start = [353, 128];    % x, y in image coord [pix]
    P_end = [260, 84];       % x, y in image coord [pix]
    scatter(P_start(:,1), P_start(:,2), 'm', 'filled')
    scatter(P_end(:,1), P_end(:,2), 'm', 'filled')
    text(P_start(:,1) + text_offs, P_start(:,2) + text_offs, 'start', 'Color', 'm')
    text(P_end(:,1) + text_offs, P_end(:,2) + text_offs, 'end', 'Color', 'm')
    
else
    error('Specify a valid method to enter channel start/end points (under "user specifiable variables" in config.m). ')
end


%% map channel geometries and extract profiles
%  and visualize on overview figure

f = figure(1); 
hold off
imagesc(DEM)
axis image
colormap gray
title('mapped channel geometries (overview)')
xlabel('x coord [pix]')
ylabel('y coord [pix]')
hold on

x_cent = cell(no_channels, 1);  
y_cent = cell(no_channels, 1); 
channel_length = cell(no_channels, 1); 
profiles = cell(no_channels, 1); 
x_prof = cell(no_channels, 1); 
y_prof = cell(no_channels, 1); 
edge_idx = cell(no_channels, 1); 
edge_coord = cell(no_channels, 1); 
edge_elev = cell(no_channels, 1); 
fchannel = cell(no_channels, 1); 

% to store whether we successfully found the channel centerline: 
channel_status = zeros(no_channels, 1); 

% loop over channels
for c = 1:no_channels
    disp(append("Start mapping channel geometry of ", channel_label(c), ". ")) 
    
    % find centerline
    [x_cent{c}, y_cent{c}, cent_length] = find_centerline(P_start(c,:), P_end(c,:), DEM, R, ... 
                                                'search_step',      search_step, ... 
                                                'search_angle',     cent_search_angle, ... 
                                                'no_samp_pts',      no_cent_samp_pts, ... 
                                                'min_diff_thr',     crack_thr, ... 
                                                'window',           window_cent); 
    channel_length{c} = sum(cent_length);       % in [pix]
    channel_length{c} = channel_length{c}*res;  % in [m]

    if isnan(channel_length{c}) ~= 1
        % found channel end
        channel_status(c) = 1;
    end

    % find cross sectional profiles
    [profiles{c}, x_prof{c}, y_prof{c}] = find_profiles(x_cent{c}, y_cent{c}, DEM, R, ...
                                                'samp_step',        prof_samp_step, ...
                                                'no_samp_pts',      no_prof_samp_pts); 
    no_profiles = size(profiles, 2); 

    % find channel edges/outlines
    window_edge = ceil(window_edge/res);        % from m to [pix]
    [edge_idx{c}, edge_coord{c}, edge_elev{c}] = find_edges(profiles{c}, x_prof{c}, y_prof{c}, prof_samp_step, slope_thr, window_edge); 
    
    % vizualise
    % centerlines
    scatter(x_cent{c}, y_cent{c}, 15, 'r', 'filled')
    plot(x_cent{c}, y_cent{c}, 'r')
    % outlines
    scatter(edge_coord{c}(:,1), edge_coord{c}(:,2), 15, 'g', 'filled')
    scatter(edge_coord{c}(:,3), edge_coord{c}(:,4), 15, 'g', 'filled')
    plot(edge_coord{c}(:,1), edge_coord{c}(:,2), 'g')
    plot(edge_coord{c}(:,3), edge_coord{c}(:,4), 'g')
    % profile transects
    scatter(x_prof{c}(:), y_prof{c}(:), 2, 'w')
    % annotation
    text(P_start(c,1) + text_offs, P_start(c,2) + text_offs, channel_label(c), 'Color', 'm')
    pause(0.01) % just to force matlab to plot
    
    % label for shapefile
    fchannel{c} = channel_label(c); 
    
end

disp(append("Finished mapping. End point reached for ", string(sum(channel_status)), "/", string(length(channel_status)), " channels. "))


%% write to files

if save_figs == 1 || save_shps == 1
    disp("Writing figure- and shapefiles... ")

    % print overview figure to file
    if save_figs == 1
        %f.WindowState = 'maximized'; % make figure fullscreen before saving
        fn = append(fig_dir, file_prefix, 'mapped_channels'); 
        print(fn, figs_filetype, figs_resolution)
        %f.WindowState = 'normal'; 
    end
    
    % write mapped geometries to shapefile
    if save_shps == 1
    
        % all centerlines in a single file
        fn = append(shp_dir, file_prefix, 'all_centerlines'); 
        lines_to_shp(x_cent, y_cent, R, fn, 'channel_label', fchannel); 
        
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

    disp(append("Done writing files. Check '", append(results_dir, proj_subdir), "' for output. "))
end


%% extended figures
%  cross sectional profiles, metrics along channel length etc.

if ext_figs == 1
disp("Creating and possibly saving extended figures. Sit tight. ")
    
    for c = 1:no_channels
        
        % cross sectional profiles
        no_profiles = size(profiles{c}, 2); 
        
        % for plotting profiles with [m] on x-axis
        prof_dist_vector = [1:size(profiles{c}, 1)]*prof_samp_step; 
        prof_dist_vector = prof_dist_vector - mean(prof_dist_vector); 

        % full cross sectional profiles using absolute elevation
        figure
        plot(prof_dist_vector, profiles{c}, 'LineWidth', 3)
        xlabel('distance from profile center [m]')
        ylabel('elevation [m]')
        title('channel cross sectional profiles (full, abs. heights)')
        % color gradient
        cmap = parula(no_profiles); 
        set(gca(), 'ColorOrder', cmap)
        hcb = colorbar; 
        title(hcb, 'norm. dist. along channel [-]')

        if save_figs == 1
            fn = append(fig_dir, file_prefix, channel_label(c), '_full_profiles_elev'); 
            print(fn, figs_filetype, figs_resolution)
        end

        % between channel edges using relative elevation (depth w.r.t. left channel edge)
        figure
        hold on
        for i = 1:no_profiles
            prof = profiles{c}(:,i); 

            % replace values outside of channel edges to NaN
            prof(1:edge_idx{c}(i,2)) = NaN;
            prof(edge_idx{c}(i,1):end) = NaN; 

            % from absolute height to depth below left channel edge
            prof = prof - edge_elev{c}(i,1); 

            plot(prof_dist_vector, prof, 'LineWidth', 3)
        end
        xlabel('distance from profile center [m]')
        ylabel('depth [m]')
        title('channel cross sectional profiles (depth below left edge)')
        cmap = parula(no_profiles); 
        set(gca(), 'ColorOrder', cmap)
        hcb = colorbar; 
        title(hcb, 'norm. dist. along profile [-]')

        if save_figs == 1
            fn = append(fig_dir, file_prefix, channel_label(c), '_lim_profiles_depth'); 
            print(fn, figs_filetype, figs_resolution)
        end
        
        
        % plotting some metrics along channel length
        norm_dist_vector = [1:no_profiles]./no_profiles; 
        
        figure(10)
        hold on
        plot(norm_dist_vector*channel_length{c}/1000, mean(profiles{c}))
        
        trough_elev = min(profiles{c}); % TO DO: min of profile not necessarily trough! could also be crack
        trough_depth = min(profiles{c})-edge_elev{c}(:,1)'; 
        channel_width = (edge_idx{c}(:,2)-edge_idx{c}(:,1))*prof_samp_step; % [m]

        figure(11)
        hold on
        plot(norm_dist_vector*channel_length{c}/1000, trough_depth)
        
        figure(12)
        hold on
        plot(norm_dist_vector*channel_length{c}/1000, channel_width/1000)  
    end
    
    figure(10)
    ylabel('elevation [m]')
    xlabel('distance along channel [km]')
    title('mean profile elevation (full) vs. distance along channel')
    legend(channel_label)
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'elev_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end
    
    
    figure(11)
    ylabel('depth w.r.t. left channel edge [m]')
    xlabel('distance along channel [km]')
    title('channel depth vs. distance along channel')
    legend(channel_label)
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'depth_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end
    
    
    figure(12)
    ylabel('channel width [km]')
    xlabel('distance along channel [km]')
    title('channel width vs. distance along channel')
    legend(channel_label)
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'width_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end
    
end
    
disp("Done! ")
