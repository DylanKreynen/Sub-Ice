clear all
close all
clc

addpath("./functions")

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
%  June - July 2024
% 
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% user specifiable variables
%  (update as required)

% path to DEM (should be GeoTIFF)
path_to_DEM = '.\input\venable.tif'; 

% output behaviour: 
results_dir = '.\output\';
proj_subdir = 'venable\'; 
fig_subdir = 'fig\'; 
shp_subdir = 'shp\'; 
file_prefix = 'default_'; % (optional)
% output path will be constructed as follows: 
% results_dir\proj_subdir\fig- or shp_subdir\file_prefix_....ext

save_figs = 1;              % print figures to disk Y/N
figs_filetype = '-dpng';    % for use with "print()"
figs_resolution = '-r500';  % for use with "print()"
ext_figs = 1;               % plot (and print) extended figures Y/N
save_shps = 1;              % save output as shapefiles Y/N

% select method to specify channel start/end points
start_end_method = 2;
% 1 = click on start/end points
% 2 = read from shapefile
% 3 = manually enter in script
path_to_start_end_shp = '.\input\venable_start_end_3.shp'; 
% ^ only needed when start_end_method is set to "2" (read from shapefile)
% important: only works if shapefile has same map projection as DEM! 

% centerline search parameters
search_step = 800;                % distance to step away from last known centerline point to construct search profile [m]
no_cent_samp_pts = 20;            % number of sampling points on search profile [-]
cent_samp_step = 100;             % distance between sampling points on search profile [m]
max_no_cent_pts = 100;            % when to stop looking for centerline end point [-]

% channel cross sectional profile parameters
prof_samp_step = 100;             % sdistance between sampling points on profile [m]
no_prof_samp_pts = 50;            % number of sampling points on profile [-]
% note: profile length ~ no_sampling_points*prof_samp_step

% slope threshold for identifying channel edge
slope_thr = 0.25;                 % [deg]


%% read DEM from GeoTIFF
%  and some basic manipulation (update as required)

[DEM, R] = readgeoraster(path_to_DEM);
res = R.CellExtentInWorldX;     % resolution of DEM [m]

% remove no data (-999), aid visualization: 
DEM(DEM<-10) = NaN; 
DEM(DEM>50) = 50; 


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

if start_end_method == 1
    % click on start/end points
    title('please specify channel start and end points, respectively:')
    P_start = ginput(1);
    scatter(P_start(1), P_start(2), 'm', 'filled')
    P_end = ginput(1);
    scatter(P_end(1), P_end(2), 'm', 'filled')
    
elseif start_end_method == 2
    % read start/end points from shapefile
    % only works if shapefile has same map projection as DEM!
    S = shaperead(path_to_start_end_shp); 
    [x_startend, y_startend] =  worldToIntrinsic(R, vertcat(S.X), vertcat(S.Y));
    P_start = [x_startend(1) y_startend(1)]; 
    P_end = [x_startend(2) y_startend(2)]; 
    scatter(P_start(1), P_start(2), 'm', 'filled')
    scatter(P_end(1), P_end(2), 'm', 'filled')

elseif start_end_method == 3
    % enter start/end points manually below
    P_start = [353, 128];    % x, y in image coord [pix]
    P_end = [260, 84];       % x, y in image coord [pix]
    scatter(P_start(1), P_start(2), 'm', 'filled')
    scatter(P_end(1), P_end(2), 'm', 'filled')
else
    error('Specify a valid method to enter channel start/end points (under "user specifiable variables"). ')
end

text_offs = 3;  % [pix]
text(P_start(1) + text_offs, P_start(2) + text_offs, 'channel start', 'Color', 'm')
text(P_end(1) + text_offs, P_end(2) + text_offs, 'channel end', 'Color', 'm')


%% find channel centerline

% basic idea: 
% - find direction (based on start/end points or previous centerline points)
% - step one search step in that direction
% - construct perpendicular sampling vector at that point
% - sample DEM (interpolate) to find elevation profile
% - new centerline location is where profile reaches min
% - repeat until channel end: step in upd dir, sample, find min
% >> see "find_centerline" function

[x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, search_step, cent_samp_step, no_cent_samp_pts, max_no_cent_pts);

% % check "anticrack" version of "find_centerline" for a crack threshold
% % > if new centerline point is significantly lower than last known
% % centerline point, we likely found a crack > select next local min instead
% crack_thr = 10; 
% [x_cent, y_cent, cent_length] = find_centerline_anticrack(P_start, P_end, DEM, R, search_step, cent_samp_step, no_cent_samp_pts, max_no_cent_pts, crack_thr);

channel_length = sum(cent_length);      % in [pix]
channel_length = channel_length*res;    % in [m]

figure(1)
hold on
scatter(x_cent, y_cent, 15, 'r', 'filled')
plot(x_cent, y_cent, 'r')

% write centerline to shapefile
if save_shps == 1
    fn = append(shp_dir, file_prefix, 'centerline'); 
    lines_to_shp(x_cent, y_cent, R, 'centnumber', 1, fn); 
end


%% find cross sectional profiles

% basic idea: 
% - find midpoint between two centerline locations 
% - find perpendicular direction to centerline
% - create sampling vector in this direction
% - sample DEM (interpolate) to find elevation profile
% - plot/store profiles
% >> see "find_profiles" function

[profiles, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, prof_samp_step, no_prof_samp_pts); 
no_profiles = size(profiles, 2); 

if save_shps == 1
    fn = append(shp_dir, file_prefix, 'profiles'); 
    lines_to_shp(x_prof, y_prof, R, 'profnumber', 1:no_profiles, fn); 
end


%% find profiles' channel edges
% now based on slope threshold, reached after reaching max slope
% >> see "find edges" function

% This scheme needs improving! Karen's idea: smooth the data using a mean
% filter, before finding the max negative curvature (see e-mail 2024-07-02)
% update "find_edges" function as required

[edge_idx, edge_coord, edge_elev] = find_edges(profiles, x_prof, y_prof, prof_samp_step, slope_thr); 

figure(1)
scatter(edge_coord(:,1), edge_coord(:,2), 15, 'g', 'filled')    % left channel edge
scatter(edge_coord(:,3), edge_coord(:,4), 15, 'g', 'filled')    % right channel edge
plot(edge_coord(:,1), edge_coord(:,2), 'g')
plot(edge_coord(:,3), edge_coord(:,4), 'g')
scatter(x_prof(:), y_prof(:), 2, 'w')

if save_figs == 1
    fn = append(fig_dir, file_prefix, 'mapped_channel'); 
    print(fn, figs_filetype, figs_resolution)
end

if save_shps == 1 
    edge_x = [edge_coord(:,1) edge_coord(:,3)];
    edge_y = [edge_coord(:,2) edge_coord(:,4)];
    fn = append(shp_dir, file_prefix, 'edges');
    lines_to_shp(edge_x, edge_y, R, 'side', ["left", "right"], fn);     
end


%% extended figures
%  cross sectional profiles, metrics along channel length etc.

if ext_figs == 1

    %%% cross sectional profile vizualisation %%%

    % for plotting profiles with [m] on x-axis
    prof_dist_vector = [1:size(profiles, 1)]*prof_samp_step; 
    prof_dist_vector = prof_dist_vector - mean(prof_dist_vector); 

    % full cross sectional profiles using absolute elevation
    figure(2)
    plot(prof_dist_vector, profiles, 'LineWidth', 3)
    xlabel('distance from profile center [m]')
    ylabel('elevation [m]')
    title('channel cross sectional profiles (full, abs. heights)')
    % color gradient
    cmap = parula(size(profiles,2)); 
    set(gca(), 'ColorOrder', cmap)
    hcb = colorbar; 
    title(hcb, 'norm. dist. along channel [-]')
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'full_profiles_elev'); 
        print(fn, figs_filetype, figs_resolution)
    end

    % between channel edges using relative elevation (depth w.r.t. left channel edge)
    figure(3)
    hold on
    for i = 1:no_profiles
        prof = profiles(:,i); 

        % replace values outside of channel edges to NaN
        prof(1:edge_idx(i,1)) = NaN; 
        prof(edge_idx(i,2):end) = NaN;

        % from absolute height to depth below left channel edge
        prof = prof - edge_elev(i,1); 

        plot(prof_dist_vector, prof, 'LineWidth', 3)
    end
    xlabel('distance from profile center [m]')
    ylabel('depth [m]')
    title('channel cross sectional profiles (depth below left edge)')
    cmap = parula(size(profiles,2)); 
    set(gca(), 'ColorOrder', cmap)
    hcb = colorbar; 
    title(hcb, 'norm. dist. along profile [-]')
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'lim_profiles_depth'); 
        print(fn, figs_filetype, figs_resolution)
    end


    %%% plotting metrics along channel length %%%

    % normalized distance along channel for every profile
    norm_dist_vector = [1:no_profiles]./no_profiles; 

    figure(4)
    hold on
    plot(norm_dist_vector*channel_length/1000, mean(profiles))
    plot(norm_dist_vector*channel_length/1000, min(profiles))
    ylabel('elevation [m]')
    xlabel('distance along channel [km]')
    legend('mean profile elevation (full profiles)', 'min. profile elevation ("centerline")', 'Location', 'southwest')
    title('mean and min. elevation vs. distance along channel')
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'elev_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end

    trough_elev = min(profiles); 
    trough_depth = min(profiles)-edge_elev(:,1)'; 
    channel_width = (edge_idx(:,2)-edge_idx(:,1))*prof_samp_step; % [m]

    figure(5)
    plot(norm_dist_vector*channel_length/1000, trough_depth)
    ylabel('depth w.r.t. left channel edge [m]')
    xlabel('distance along channel [km]')
    title('channel depth vs. distance along channel')
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'depth_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end

    figure(6)
    plot(norm_dist_vector*channel_length/1000, channel_width/1000)
    ylabel('channel width [km]')
    xlabel('distance along channel [km]')
    title('channel width vs. distance along channel')
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'width_vs_distance_along_channel'); 
        print(fn, figs_filetype, figs_resolution)
    end


    %%% to check channel edges for a single profile %%%

    i = 15; 
    prof = profiles(:,i); 
    l_idx = edge_idx(i,1);
    r_idx = edge_idx(i,2);
    slope = abs(rad2deg(gradient(prof))); 

    figure
    hold on
    plot(prof)
    plot(slope)
    scatter(l_idx, prof(l_idx), 'g')
    scatter(r_idx, prof(r_idx), 'r')
    scatter(l_idx, slope(l_idx), 'g')
    scatter(r_idx, slope(r_idx), 'r')
    legend('prof. elev. [m]', 'abs. slope [deg]', 'left edge [-]', 'right edge [-]')
    title('example cross sectional profile incl. identified channel edges')
    xlabel('along profile distance in sampling steps [-]')
    ylabel('elevation [m]')
    
    if save_figs == 1
        fn = append(fig_dir, file_prefix, 'example_profile_incl_edges'); 
        print(fn, figs_filetype, figs_resolution)
    end

end

