function [x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, varargin) 
%[x_cent, y_cent, cent_length] = find_centerline_anticrack(P_start, P_end, DEM, R, search_step, search_angle, no_samp_pts, max_no_cent_pts, min_diff_thr, window)
%Returns the coordinates of the channel centerline. 
% basic idea: 
% - find direction (based on start/end points or previous centerline points)
% - step one search step in that direction
% - construct circular arc series of sampling points
% - sample DEM (interpolate) to find elevation profile on circular arc
% - find local min and select appropriate one as new centerline point
% - repeat until channel end: step in upd dir, sample, find min
% - return centerline coordinates and section lengths
%
% required input:
% P_start = vector containing x and y img coordinates of start point [pix]
% P_end = vector containing x and y ing coordinates of end point [pix]
% DEM = elevation data array [m] (use readgeoraster to read a geotiff)
% R = spatial referencing information for the array [-]
% 
% optional input: 
% search_step = distance to step away from P_start to construct search profile [m] (default: 1000m)
% search_angle = angle of view within to look for centerline [deg] (default: 60deg)
% no_samp_pts = number of sampling points on search profile [-] (default: 10)
% max_no_cent_pts = when to stop looking for centerline end point [-] (default: 50m)
% min_diff_thr = if new centerline location has an elevation value that's 
%                significantly different from the previously known 
%                centerline location (determined by threshold in [m]), we
%                likely found a crack and select the next local min instead (default: 100m)
% window = window size for search profile smoothing [pix], set to 0 for no smoothing (default)
%
% output:
% x_cent = vector with x coordinates of channel centerline [pix]
% y_cent = vector with y coordinates of channel centerline [pix]
% cent_length = vector with centerline section lengths [pix]
% 
% (c) Dylan Kreynen
% University of Oslo
% June - Oct 2024


%% inputParser

% default parameter values
default_search_step = 1000; 
default_search_angle = 60; 
default_no_samp_pts =  10; 
default_max_no_cent_pts = 50; 
default_min_diff_thr = 100; 
default_window = 0; 

% parse input arguments
p = inputParser; 
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validMapCellsRef = @(x) class(x) == "map.rasterref.MapCellsReference"; 
validPStartEnd = @(x) (size(x, 1)==2 | size(x, 2)==2) && isnumeric(x(1)) && isscalar(x(1)); 
addRequired(p, 'P_start', validPStartEnd)
addRequired(p, 'P_end', validPStartEnd)
addRequired(p, 'DEM')
addRequired(p, 'R', validMapCellsRef)
addOptional(p, 'search_step', default_search_step, validScalarPosNum)
addOptional(p, 'search_angle', default_search_angle, validScalarPosNum)
addOptional(p, 'no_samp_pts', default_no_samp_pts, validScalarPosNum)
addOptional(p, 'max_no_cent_pts', default_max_no_cent_pts, validScalarPosNum)
addOptional(p, 'min_diff_thr', default_min_diff_thr, validScalarPosNum)
addOptional(p, 'window', default_window, validScalarPosNum)
parse(p, P_start, P_end, DEM, R, varargin{:}); 

search_step = p.Results.search_step; 
search_angle = p.Results.search_angle; 
no_samp_pts = p.Results.no_samp_pts; 
max_no_cent_pts = p.Results.max_no_cent_pts; 
min_diff_thr = p.Results.min_diff_thr; 
window = p.Results.window; 


%% actual function

res = R.CellExtentInWorldX;        % resolution of DEM [m/pix]
DEM_int = griddedInterpolant(DEM); % interpolant for our DEM

% stop criterion: 
stop_dist = search_step; 
stop_dist = stop_dist/res;  % now in [pix]
% when we are less than a search step away from the end point, we have
% reached our destination/centerline is considered complete

% preallocate for speed: 
x_cent = NaN(max_no_cent_pts, 1); 
y_cent = NaN(max_no_cent_pts, 1);
cent_length = NaN(max_no_cent_pts-1, 1); 

% start searching; find first point after "start": 
% construct sampling vector, one search step in the direction of channel end 
[x_samp, y_samp] = centerline_query_pts(P_start, P_end, res, search_step, search_angle, no_samp_pts, 1); 
% sample DEM at this sampling vector to get a profile
prof = DEM_int(y_samp, x_samp); 
% smooth profile using mean filter
if window ~= 0
    prof = smoothdata(prof, 'movmean', window); 
end
% find min of profile = new centerline location
[min_val, min_loc] = min(prof);

% add start and first point to centerline: 
x_cent(1:2) = [P_start(1); x_samp(min_loc)]; 
y_cent(1:2) = [P_start(2); y_samp(min_loc)]; 

% store section length: 
cent_length(1) = sqrt((x_cent(1)-x_cent(2))^2 + (y_cent(1)-y_cent(2))^2); 

% find the next centerline points (until close to end point)
no_pts = ceil(no_samp_pts/2);
middle_weights = [1:no_pts, no_pts+1, abs(-no_pts:-1)]; 
% ^ vector with weights - higher closer to middle of profile
i = 2; 
dist_to_end = stop_dist + 1; 
while dist_to_end > stop_dist % give condition here (distance to end point)
    i = i + 1; 
    if i > max_no_cent_pts-1
        disp("Warning: max. no. of centerline points reached, but did not reach channel end.")
        break
    end

    [x_samp, y_samp] = centerline_query_pts([x_cent(i-2) y_cent(i-2)], [x_cent(i-1) y_cent(i-1)], res, search_step, search_angle, no_samp_pts, 0);
    prof = DEM_int(y_samp, x_samp); % evaluate interpolant of 
    % smooth profile using mean filter
    if window ~= 0
        prof = smoothdata(prof, 'movmean', window); 
    end
    
    % update "anti-crack": rather than find the absolute min of the profile,
    % we find all local min. then select the appropriate one
    %  - local min which best preserves centerline direction, 
    %  - whose elev diff wrt prev centerline point does not exceed threshold
    % if no appropriate local min can be found we preserve direction from
    % previous centerline section

    % another idea: basal channels tend to be smooth, cracks sharp
    % can we use that to distinguish between the channel and a crack? 
    % e.g. make cross sectional profiles on the go and compare slopes

    % another idea: rather than step in the direction of the last
    % centerline section, we could step in the mean direction of the last x
    % number of centerline sections (less likely to get stuck in a crack?)

    % one more idea: we could sample a profile between the last known
    % centerline point and the new potential centerline point. If it's
    % indeed along the centerline we should not have big bumps in it. If
    % it's a crack, there's likely "a big hill" on the profile

    % find all local min on profile
    local_min_loc = islocalmin(prof); % logical

    while any(local_min_loc)   % while we have valid local minima
        % select local min closest to center of search profile
        local_min_loc = local_min_loc.*middle_weights; 
        [~, min_loc] = max(local_min_loc); 
        min_val_upd = prof(min_loc); 

        % check for depth wrt previous centerline point threshold
        if (min_val - min_val_upd) > min_diff_thr
            % disp(append("Avoided a CRACK! Cent. idx. ", string(i-1)))
            % remove local min from list and try again
            local_min_loc(min_loc) = 0; 
        else
            % we found an appropriate local min
            % break out of while loop
            break
        end
    end

    if ~any(local_min_loc)
        % no local min on profile that satisfies conditions
        % > stick to direction from previous section
        min_loc = no_pts+1; % center of profile
        min_val_upd = prof(min_loc); 
    end

    min_val = min_val_upd; 

    % new centerline point: 
    x_cent(i) = x_samp(min_loc); 
    y_cent(i) = y_samp(min_loc); 
    
    % store section length: 
    cent_length(i-1) = sqrt((x_cent(i)-x_cent(i-1))^2 + (y_cent(i)-y_cent(i-1))^2); 
    
    % compute distance to end point: 
    dist_to_end = sqrt((x_cent(i)-P_end(1))^2 + (y_cent(i)-P_end(2))^2); 
end

% add user specified end point to centerline: 
x_cent(i+1) = P_end(1); 
y_cent(i+1) = P_end(2); 

% strip of unnecessary NaNs: 
x_cent = x_cent(1:i+1);
y_cent = y_cent(1:i+1); 

% compute and store last section lengths: 
cent_length(i) = sqrt((x_cent(end-1)-x_cent(end))^2 + (y_cent(end-1)-y_cent(end))^2);
cent_length = cent_length(1:i); 

end