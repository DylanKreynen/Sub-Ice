function [x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, search_step, samp_step, no_samp_pts, max_no_cent_pts)
%[x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, search_step, samp_step, no_samp_pts, max_no_cent_pts)
%Returns the coordinates of the channel centerline. 
% basic idea: 
% - find direction (based on start/end points or previous centerline points)
% - step one search step in that direction
% - construct perpendicular sampling vector at that point
% - sample DEM (interpolate) to find elevation profile
% - new centerline location is where profile reaches min
% - repeat until channel end: step in upd dir, sample, find min
% - return centerline coordinates and section lengths
%
% input: 
% P_start = vector containing x and y img coordinates of start point [pix]
% P_end = vector containing x and y img coordinates of end point [pix]
% DEM = elevation data array [m] (use readgeoraster to read a geotiff)
% R = spatial referencing information for the array [-]
% search_step = distance to step away from P_start to construct search profile [m]
% samp_step = distance between sampling points on search profile [m]
% no_samp_pts = number of sampling points on search profile [-]
% max_no_cent_pts = when to stop looking for centerline end point [-]
%
% output: 
% x_cent = vector with x coordinates of channel centerline [pix]
% y_cent = vector with y coordinates of channel centerline [pix]
% cent_lengths = vector with centerline section lengths [pix]
% 
% (c) Dylan Kreynen
% University of Oslo
% June - July 2024

x_start = P_start(1); 
y_start = P_start(2); 
x_end = P_end(1); 
y_end = P_end(2); 

res = R.CellExtentInWorldX; % resolution of DEM [m/pix]
x = 1:size(DEM, 1); 
y = 1:size(DEM, 2);
[X, Y] = meshgrid(y, x);

% stop criterion: 
stop_dist = search_step; 
stop_dist = stop_dist/res;  % now in [pix]
% when we are less than a search step away from the end point, we have
% reached our destination/centerline is considered complete

% start searching; find first point after "start": 
% construct sampling vector, one search step in the direction of channel end 
[x_samp, y_samp] = centerline_query_pts(P_start, P_end, res, search_step, samp_step, no_samp_pts, 1); 
% sample DEM at this sampling vector to get a profile
profile = interp2(X, Y, DEM, x_samp, y_samp); 
% find min of profile = new centerline location
[~, min_loc] = min(profile);

% add start and first point to centerline: 
x_cent = [x_start; x_samp(min_loc)]; 
y_cent = [y_start; y_samp(min_loc)]; 

% store section length: 
cent_length(1) = sqrt((x_cent(1)-x_cent(2))^2 + (y_cent(1)-y_cent(2))^2); 

% find the next centerline points (until close to end point)
i = 2; 
dist_to_end = stop_dist + 1; 
while dist_to_end > stop_dist % give condition here (distance to end point)
    i = i + 1; 
    if i > max_no_cent_pts-1
        disp("Warning: max. no. of centerline points reached, but did not reach channel end.")
        break
    end

    x1 = x_cent(i-2); y1 = y_cent(i-2);     % second most recent centerline point coords
    x2 = x_cent(i-1); y2 = y_cent(i-1);     % most recent centerline point coords
    
    [x_samp, y_samp] = perp_search_sampling(x1, y1, x2, y2, res, search_step, samp_step, no_samp_pts, 0);
    profile = interp2(X, Y, DEM, x_samp, y_samp); 
    [~, min_loc] = min(profile);
    
    % new centerline point is wherever sampled elevation profile is lowest: 
    xn = x_samp(min_loc); x_cent(i) = xn; 
    yn = y_samp(min_loc); y_cent(i) = yn; 
    
    % store section length: 
    cent_length(i-1) = sqrt((xn-x2)^2 + (yn-y2)^2); 
    
    % compute distance to end point: 
    dist_to_end = sqrt((xn-x_end)^2 + (yn-y_end)^2); 
end

% add user specified end point to centerline: 
x_cent(end+1) = x_end; 
y_cent(end+1) = y_end; 

% compute and store last section lengths: 
cent_length(end+1) = sqrt((x_cent(end-1)-x_cent(end))^2 + (y_cent(end-1)-y_cent(end))^2);

end

