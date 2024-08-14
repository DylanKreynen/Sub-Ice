clear all
close all
clc

%% Sub-Ice: DEM-based semi-automized mapping of ice shelf basal channels
%  Algorithm that returns ice shelf basal channels' centerline, outlines
%  and cross sectional profiles; based on surface expressions in a DEM. 
%  This script maps a single channel over different DEMs, given user
%  identiefied channel start and end points. (See "map_multi_channel.m" for
%  mapping one or more channels given a single DEM.)
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
%  June - August 2024
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


%% create output directories

% figures
if save_figs == 1
    fig_dir = append(results_dir, proj_subdir, fig_subdir); 
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir)
    end
end

% shapefiles
if save_shps == 1
    shp_dir = append(results_dir, proj_subdir, shp_subdir); 
    if ~exist(shp_dir, 'dir')
        mkdir(shp_dir)
    end
end


%% list DEMs/GeoTIFFs
%  and corresponding dates

DEM_files = dir(append(path_to_DEMs, '*.tif')); 
DEM_files = {DEM_files.name}'; 
no_DEMs = length(DEM_files); 

dates = NaT(no_DEMs,1); 
for f = 1:no_DEMs
    yyyymmdd = DEM_files{f}(end-11:end-4); 
    yyyymmdd = double(string(yyyymmdd)); 
    dates(f) = datetime(yyyymmdd, 'ConvertFrom', 'yyyymmdd'); 
end

% if files not sorted chronologically, sort
if ~issorted(dates)
    [dates, idx] = sort(dates); 
    DEM_files = DEM_files(idx); 
end


%% specify channel centerline start/end points
%  have user click on channel start and end points through GUI, read start
%  and end points from shapefile or enter img coordinates manually 

text_offs = 15;  % [pix] (for plotting only)

if start_end_method == 1
    % click start/end points
    P_start = NaN(no_DEMs, 2); 
    P_end = NaN(no_DEMs, 2); 
    
    for t = 1:no_DEMs    
        DEM = readgeoraster(append(path_to_DEMs, DEM_files{t})); 
        DEM(DEM<-10) = NaN;  % remove no data (-9999), 
        DEM(DEM>50) = 50;    % and aid vizualisation
    
        figure(1)
        hold off
        imagesc(DEM)
        hold on
        axis image
        colormap gray
        xlabel('x coord [pix]')
        ylabel('y coord [pix]')
        
        title(append(string(dates(t)), " - please left click on channel start and end points, respectively."))
        disp(append(string(dates(t)), " - please left click on channel start and end points, respectively."))
        [P_start(t,1), P_start(t,2)] = ginput(1); 
        scatter(P_start(t,1), P_start(t,2), 'm', 'filled')
        text(P_start(t,1) + text_offs, P_start(t,2) + text_offs, 'channel start', 'Color', 'm')  
        [P_end(t,1), P_end(t,2)] = ginput(1); 
        text(P_end(t,1) + text_offs, P_end(t,2) + text_offs, 'channel end', 'Color', 'm')
        scatter(P_end(t,1), P_end(t,2), 'm', 'filled')
        pause(0.5)
    end
    pause(0.5)
    close figure 1

elseif start_end_method == 2
    % read start/end points from shapefile
    % only works if shapefile has same map projection as first DEM!
    [~, R] = readgeoraster(append(path_to_DEMs, DEM_files{1})); 
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

elseif start_end_method == 3
    % enter start/end points manually below
    P_start = [353, 128];    % x, y in image coord [pix]
    P_end = [260, 84];       % x, y in image coord [pix]
    
else
    error('Specify a valid method to enter channel start/end points (under "user specifiable variables"). ')
end

% check if we have start/end points for all DEMs
no_start_end_pairs = size(P_start(~isnan(P_start)), 1)/2; 
if no_start_end_pairs ~= no_DEMs
    error('Number of start/end point pairs does not match number of DEMs. ')
end


%% map channel geometries and extract profiles
%  and visualize on overview figure

x_cent = cell(no_DEMs, 1);  
y_cent = cell(no_DEMs, 1); 
channel_length = cell(no_DEMs, 1); 
profiles = cell(no_DEMs, 1); 
x_prof = cell(no_DEMs, 1); 
y_prof = cell(no_DEMs, 1); 
edge_idx = cell(no_DEMs, 1); 
edge_coord = cell(no_DEMs, 1); 
edge_elev = cell(no_DEMs, 1); 
fchannel = cell(no_DEMs, 1); 

% loop over timesteps
for t = 1:no_DEMs
    disp(append("Start mapping channel geometry on ", string(dates(t)), ". "))

    % read DEM from GeoTIFF
    [DEM, R] = readgeoraster(append(path_to_DEMs, DEM_files{t}));
    res = R.CellExtentInWorldX;     % resolution of DEM [m]
    DEM(DEM<-10) = NaN;  % remove no data (-999), 
    DEM(DEM>50) = 50;    % and aid vizualisation (update as required)
    
    % find channel centerline and centerline length
    window_cent = ceil(window_cent/res);        % from m to [pix]
    [x_cent{t}, y_cent{t}, cent_length] = find_centerline(P_start(t,:), P_end(t,:), DEM, R, ... 
                                                'search_step',      search_step, ... 
                                                'samp_step',        cent_samp_step, ... 
                                                'no_samp_pts',      no_cent_samp_pts, ... 
                                                'max_no_cent_pts',  max_no_cent_pts, ... 
                                                'min_diff_thr',     crack_thr, ... 
                                                'window',           window_cent); 

    channel_length{t} = sum(cent_length);       % in [pix]
    channel_length{t} = channel_length{t}*res;  % in [m]

    % find cross sectional profiles
    [profiles{t}, x_prof{t}, y_prof{t}] = find_profiles(x_cent{t}, y_cent{t}, DEM, R, ...
                                                'samp_step',        prof_samp_step, ...
                                                'no_samp_pts',      no_prof_samp_pts); 
    no_profiles = size(profiles, 2); 

    % find channel edges/outlines
    window_edge = ceil(window_edge/res);        % from m to [pix]
    [edge_idx{t}, edge_coord{t}, edge_elev{t}] = find_edges(profiles{t}, x_prof{t}, y_prof{t}, prof_samp_step, slope_thr, window_edge); 
    
    % label for shapefile
    fchannel{t} = string(dates(t)); 

    % vizualise
    figure(t)
    imagesc(DEM)
    axis image
    colormap gray
    title(append('mapped channel geometry on ', fchannel{t}))
    xlabel('x coord [pix]')
    ylabel('y coord [pix]')
    hold on
    % centerlines
    scatter(x_cent{t}, y_cent{t}, 15, 'r', 'filled')
    plot(x_cent{t}, y_cent{t}, 'r')
    % outlines
    scatter(edge_coord{t}(:,1), edge_coord{t}(:,2), 15, 'g', 'filled')
    scatter(edge_coord{t}(:,3), edge_coord{t}(:,4), 15, 'g', 'filled')
    plot(edge_coord{t}(:,1), edge_coord{t}(:,2), 'g')
    plot(edge_coord{t}(:,3), edge_coord{t}(:,4), 'g')
    % profile transects
    scatter(x_prof{t}(:), y_prof{t}(:), 2, 'w')

    pause(0.01) % just to force matlab to plot

    % print fig to file
    if save_figs == 1
        %f.WindowState = 'maximized'; % make figure fullscreen before saving
        fn = append(fig_dir, file_prefix, fchannel{t}); 
        print(fn, figs_filetype, figs_resolution)
        %f.WindowState = 'normal'; 
    end
end

disp("Finished mapping all timesteps! ")


%% write shapefiles

% write mapped geometries to shapefile
if save_shps == 1

    disp("Writing shapefiles... ")

    % all centerlines in a single file
    fn = append(shp_dir, file_prefix, 'all_centerlines'); 
    lines_to_shp(x_cent, y_cent, R, fn, 'channel_label', fchannel); 
    
    % centerlines and outlines in one file per channel
    for t = 1:no_DEMs
        outlines_x = cell(3, 1); 
        outlines_y = cell(3, 1); 
        outlines_x{1} = x_cent{t};            % centerline
        outlines_y{1} = y_cent{t}; 
        outlines_x{2} = edge_coord{t}(:,1);   % left edge
        outlines_y{2} = edge_coord{t}(:,2); 
        outlines_x{3} = edge_coord{t}(:,3);   % right edge
        outlines_y{3} = edge_coord{t}(:,4); 
        fline = {"centerline", "left_edge", "right_edge"}; 
        fn = append(shp_dir, file_prefix, fchannel{t}, "_outlines"); 
        lines_to_shp(outlines_x, outlines_y, R, fn, 'line_type', fline);
    end
    
    % all profile transects in one file per channel
    for t = 1:no_DEMs
        no_profiles = size(x_prof{t}, 2); 
        fprof = 1:no_profiles;  
        fn = append(shp_dir, file_prefix, fchannel{t}, "_profiles"); 
        lines_to_shp(x_prof{t}, y_prof{t}, R, fn, 'prof_no', fprof);
    end 
end

if save_figs == 1 | save_shps == 1
    disp(append("Done writing files. Check '", append(results_dir, proj_subdir), "' for output. "))
end


%% extended figures
%  cross sectional profiles, metrics along channel length etc.

if ext_figs == 1
disp("Creating and possibly saving extended figures. Sit tight. ")
    
    for t = 1:no_DEMs
        
        % cross sectional profiles
        no_profiles = size(profiles{t}, 2); 
        
        % for plotting profiles with [m] on x-axis
        prof_dist_vector = [1:size(profiles{t}, 1)]*prof_samp_step; 
        prof_dist_vector = prof_dist_vector - mean(prof_dist_vector); 

        % full cross sectional profiles using absolute elevation
        figure
        plot(prof_dist_vector, profiles{t}, 'LineWidth', 3)
        xlabel('distance from profile center [m]')
        ylabel('elevation [m]')
        title('channel cross sectional profiles (full, abs. heights)')
        % color gradient
        cmap = parula(no_profiles); 
        set(gca(), 'ColorOrder', cmap)
        hcb = colorbar; 
        title(hcb, 'norm. dist. along channel [-]')

        if save_figs == 1
            fn = append(fig_dir, file_prefix, fchannel{t}, '_full_profiles_elev'); 
            print(fn, figs_filetype, figs_resolution)
        end

        % between channel edges using relative elevation (depth w.r.t. left channel edge)
        figure
        hold on
        for i = 1:no_profiles
            prof = profiles{t}(:,i); 

            % replace values outside of channel edges to NaN
            prof(1:edge_idx{t}(i,2)) = NaN;
            prof(edge_idx{t}(i,1):end) = NaN; 

            % from absolute height to depth below left channel edge
            prof = prof - edge_elev{t}(i,1); 

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
            fn = append(fig_dir, file_prefix, fchannel{t}, '_lim_profiles_depth'); 
            print(fn, figs_filetype, figs_resolution)
        end
        
        
        % plotting some metrics along channel length
        norm_dist_vector = [1:no_profiles]./no_profiles; 
        
        figure(10)
        hold on
        plot(norm_dist_vector*channel_length{t}/1000, mean(profiles{t}))
        
        trough_elev = min(profiles{t}); 
        trough_depth = min(profiles{t})-edge_elev{t}(:,1)'; 
        channel_width = (edge_idx{t}(:,2)-edge_idx{t}(:,1))*prof_samp_step; % [m]

        figure(11)
        hold on
        plot(norm_dist_vector*channel_length{t}/1000, trough_depth)
        
        figure(12)
        hold on
        plot(norm_dist_vector*channel_length{t}/1000, channel_width/1000)  
    end
    
    figure(10)
    ylabel('elevation [m]')
    xlabel('distance along channel [km]')
    title('mean profile elevation (full) vs. distance along channel')
    legend(fchannel)
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'elev_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end
    
    
    figure(11)
    ylabel('depth w.r.t. left channel edge [m]')
    xlabel('distance along channel [km]')
    title('channel depth vs. distance along channel')
    legend(fchannel)
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'depth_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end
    
    
    figure(12)
    ylabel('channel width [km]')
    xlabel('distance along channel [km]')
    title('channel width vs. distance along channel')
    legend(fchannel)
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'width_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end
    
end
    
disp("Done! ")
