clear all
% close all
clc

addpath("./functions")

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
%  June - July 2024
% 
%  originally a project at the Int. Summer School in Glaciology
%  project team members: Marcelo Santis & Dylan Kreynen
%  advisor: Karen Alley (University of Manitoba)
%  McCarthy (AK), June 2024


%% user specifiable variables
%  (update as required)

% path to DEM (should be GeoTIFF)
path_to_DEM = '..\REMA\venable_10m.tif'; 
path_to_DEM = '.\input\venable.tif'; 

% output behaviour: 
results_dir = '.\output\';
proj_subdir = 'venable\'; 
fig_subdir = 'fig\'; 
shp_subdir = 'shp\'; 
file_prefix = 'default_'; % (optional)
% output path will be constructed as follows: 
% results_dir\proj_subdir\fig_subdir\file_prefix_....ext

save_figs = 0;              % print figures to disk Y/N
figs_filetype = '-dpng';    % for use with "print()"
figs_resolution = '-r500';  % for use with "print()"
ext_figs = 0;               % plot (and print) extended figures Y/N
save_shps = 0;              % save output as shapefiles Y/N

% select method to specify channel start/end points
start_end_method = 2;
% 1 = click on start/end points
% 2 = read from shapefile
% 3 = manually enter in script
path_to_start_end_shp = '.\input\venable_start_end_4.shp'; 
% ^ only needed when start_end_method is set to "2" (read from shapefile)
% important: only works if shapefile has same map projection as DEM! 

% centerline search parameters
search_step = 1000;               % distance to step away from last known centerline point to construct search profile [m]
no_cent_samp_pts = 30;            % number of sampling points on search profile [-]
cent_samp_step = 50;              % distance between sampling points on search profile [m]
max_no_cent_pts = 50;             % when to stop looking for centerline end point [-]
crack_thr = 8;                    % if new centerline point's depth w.r.t. last known point is greater than threshold, pick next best point instead [m]
window_cent = 0;                  % window size for search profile smoothing [m] (will be rounded up to [pix], set to 0 for no smoothing)

% two strategies:   1. less dense sampling (e.g. 50m step), no smoothing
%                   2. dense sampling + smoothing (e.g. 10m step, 50m smoothing)



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


%% read DEM from GeoTIFF
%  and some basic manipulation

[DEM, R] = readgeoraster(path_to_DEM);
res = R.CellExtentInWorldX;     % resolution of DEM [m]

% remove no data values, aid visualization
% (update as required)
DEM(DEM<-15) = NaN; % no data: -9999
DEM(DEM>50) = 50; 


%% specify channel centerline start/end points
%  have user click on channel start and end points through GUI, read start
%  and end points from shapefile or enter img coordinates manually

text_offs = 3;  % [pix]

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

no_channels = size(P_start, 1); 
disp(append("Found ", string(no_channels), " channels' start and end points in shapefile. "))

% channel label - feel free to adjust: 
channel_label = append("channel-", string(1:no_channels)); 
% (should be string)


%% map channel geometries and extract profiles
%  and visualize on overview figure

figure 
hold off
imagesc(DEM)
axis image
colormap gray
title('mapped channel centerlines (overview)')
xlabel('x coord [pix]')
ylabel('y coord [pix]')
hold on

x_cent = cell(no_channels, 1);  
y_cent = cell(no_channels, 1); 
channel_length = cell(no_channels, 1); 

% loop over channels
tic
for c = 1:no_channels
    disp(append("Start mapping channel centerline of ", channel_label(c), ". "))
    
    % find channel centerline and centerline length
    window_cent = ceil(window_cent/res);        % from m to [pix]
    [x_cent{c}, y_cent{c}, cent_length] = find_centerline(P_start(c,:), P_end(c,:), DEM, R, ... 
                                                'search_step',      search_step, ... 
                                                'samp_step',        cent_samp_step, ... 
                                                'no_samp_pts',      no_cent_samp_pts, ... 
                                                'max_no_cent_pts',  max_no_cent_pts, ... 
                                                'min_diff_thr',     crack_thr, ... 
                                                'window',           window_cent); 
    channel_length{c} = sum(cent_length);       % in [pix]
    channel_length{c} = channel_length{c}*res;  % in [m]
    
    % vizualise
    % centerlines
    scatter(x_cent{c}, y_cent{c}, 15, 'r', 'filled')
    plot(x_cent{c}, y_cent{c}, 'r')
    % annotation
    text(P_start(c,1) + text_offs, P_start(c,2) + text_offs, channel_label(c), 'Color', 'm')
    pause(0.05) % just to force matlab to plot
    
    % label for shapefile
    fchannel{c} = channel_label(c); 
end

disp("Finished mapping all channels! ")
toc


%% write to files

if save_figs == 1 | save_shps == 1
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
    end

    disp(append("Done writing files. Check '", append(results_dir, proj_subdir), "' for output. "))
end

