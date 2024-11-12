function [x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, varargin) 
%[x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, search_step, search_angle, max_gradient, window, max_length_factor, max_recursions)
%Returns the coordinates of the channel centerline. 
% basic idea: 
% - find direction (based on start/end points or previous centerline points)
% - step one search step in that direction
% - construct circular arc series of sampling points
% - sample DEM (interpolate) to find elevation profile on circular arc
% - find local min and select appropriate one as new centerline point
% - repeat until channel end: step in upd dir, sample, find min
% - return centerline coordinates and section lengths
% note: if the channel end point is not reached, this function can try
% again recursively with updated search parameters (set max_recursions > 1)
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
% max_gradient = if newly found centerline segment exceeds this gradient, we assume we found
%                a crack instead > take next best candidate as new centerline point [%]
% window = window size for search profile smoothing [m] (will be rounded, set to 0 for no smoothing)
% max_length_factor = controls when to stop looking for channel end point [-] (default: 1.75)
%                     max. channel length = (max_length_factor)*(distance between start and end point)
% max_recursions = max. no. of attempts with updated search parameters [-] (set to 1 for no recursion) (default = 1)
%
% output:
% x_cent = vector with x coordinates of channel centerline [pix]
% y_cent = vector with y coordinates of channel centerline [pix]
% cent_length = vector with centerline section lengths [pix]
% 
% (c) Dylan Kreynen
% University of Oslo
% 2024


%% inputParser

% default parameter values
default_search_step = 1000; 
default_search_angle = 60; 
default_max_length_factor = 1.75; 
default_max_gradient = 10; 
default_window = 0; 
default_max_recursions = 1; 

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
addOptional(p, 'max_length_factor', default_max_length_factor, validScalarPosNum)
addOptional(p, 'max_gradient', default_max_gradient, validScalarPosNum)
addOptional(p, 'window', default_window, validScalarPosNum)
addOptional(p, 'max_recursions', default_max_recursions, validScalarPosNum)
parse(p, P_start, P_end, DEM, R, varargin{:}); 

search_step = p.Results.search_step; 
search_angle = p.Results.search_angle; 
max_length_factor = p.Results.max_length_factor; 
max_gradient = p.Results.max_gradient; 
window = p.Results.window; 
max_recursions = p.Results.max_recursions; 


%% actual function

res = R.CellExtentInWorldX;        % resolution of DEM [m/pix]
DEM_int = griddedInterpolant(DEM); % interpolant for our DEM

% length of circular search segment [m]
searchlen = (search_angle/360)*2*pi*search_step; 
% no. of sampling pts on profile to roughly match DEM res [-]
no_samp_pts = ceil(searchlen/res); 
% no. of search segment pts in smoothing window [-]
window = ceil(window/res);

% stop criterion: 
stop_dist = search_step; 
stop_dist = stop_dist/res;  % now in [pix]
% when we are less than a search step away from the end point, we have
% reached our destination/centerline is considered complete

% translate max_gradient to max. elevation difference
% if newly found centerline point's elevation exceeds this value, we assume 
% we found a crack instead and take the next best local min on search profile
max_diff_elev = (max_gradient/100)*search_step; 

% alternative stop criterion: max. centerline length reached
if max_length_factor < 1
    disp("Warning: max_length_factor < 1, will not be able to reach channel centerline end point. ")
end
dist_start_end = sqrt((P_start(1)-P_end(1))^2 + (P_start(2)-P_end(2))^2);   % [pix]
dist_start_end = dist_start_end*res;                                        % [m]
max_no_cent_pts = ceil((max_length_factor*dist_start_end)/search_step);     % [-] max. no of centerline points

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
    prof = smoothdata(prof, 'movmean', window, 'omitnan'); 
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
        %disp("Warning: did not reach channel end (exceeded max. centerline length). ")
        break
    end

    [x_samp, y_samp] = centerline_query_pts([x_cent(i-2) y_cent(i-2)], [x_cent(i-1) y_cent(i-1)], res, search_step, search_angle, no_samp_pts, 0);
    prof = DEM_int(y_samp, x_samp); % evaluate interpolant of 
    % smooth profile using mean filter
    if window ~= 0
        prof = smoothdata(prof, 'movmean', window);
    end
    
    % another idea: basal channels tend to be smooth, cracks sharp
    % can we use that to distinguish between the channel and a crack? 
    % e.g. make cross sectional profiles on the go and compare slopes

    % another idea: rather than step in the direction of the last
    % centerline section, we could step in the mean direction of the last x
    % number of centerline sections (less likely to get stuck in a crack?)

    % one more idea: we could sample a profile between the last known
    % centerline point and the new potential centerline point. If it's
    % indeed along the centerline we should not have big bumps in it. If
    % it's a crack, there's likely "a big hill" on the profile (or a steep one)

    % idea: if we don't end up and the channel's predefined end point, we
    % could try to find it by reversing the start and end points (i.e. look
    % for the start point, given the end point). one more shot

    % similarly: if we don't reach the channel's end point, we could go
    % back to the step where there were multiple centerline options
    % (multiple minima) and go for the next best one and try again
    % follow up: mark every centerline point with a confidence score? 

    % observation: for a decent number of channels, things go wrong at the
    % very first step because the channel is locally orientated in quite a
    % different direction than the channel end point. idea: bigger search
    % squat for first point (180deg?)

    % find all local min on profile
    local_min_loc = islocalmin(prof); % logical

    while any(local_min_loc)   % while we have valid local minima
        % select local min closest to center of search profile
        local_min_loc = local_min_loc.*middle_weights; 
        [~, min_loc] = max(local_min_loc); 
        min_val_upd = prof(min_loc); 

        % check for depth wrt previous centerline point threshold
        if (min_val - min_val_upd) > max_diff_elev
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


%% keep trying! 
it_no = 1; 
if max_recursions > 1
    % we should check whether we reached channel end, or if we need to try
    % again with slightly different search parameters
    
    channel_length = sum(cent_length); % [pix]
    if ~isnan(channel_length)
        % found channel end, exit function
        disp("Found channel centerline on first try. ")
        return
    else
        % flip start and end points and try again
        it_no = it_no + 1; 
        [x_cent, y_cent, cent_length] = find_centerline(P_end, P_start, DEM, R, ... 
                                                            'search_step',      search_step, ... 
                                                            'search_angle',     search_angle, ... 
                                                            'max_gradient',     max_gradient, ... 
                                                            'window',           window, ...
                                                            'max_length_factor',max_length_factor); 
        channel_length = sum(cent_length); % [pix]
        if ~isnan(channel_length)
            % found channel end, flip and exit function
            x_cent = flip(x_cent); 
            y_cent = flip(y_cent); 
            cent_length = flip(cent_length); 
            disp("Found channel centerline on second try (reversed). ")
            return
        end
    end
    
    % if we end up here, we still did not find channel end... 
    % start disturbing the centerline search parameters (slightly)
    while it_no < max_recursions
        % jiggle input parameters
        pret = randn(5, 1);      % preturbations, randomly drawn from normal distribution
        pret = 1 + pret*0.05;    % rescale: small preturbations only
        loop_search_step = pret(1)*search_step;
        loop_search_angle = pret(2)*search_angle;
        loop_max_gradient = pret(4)*max_gradient;
        loop_window = pret(5)*window;
    
        % try and find centerline again (start to end and end to start)
        it_no = it_no + 1; 
        [x_cent, y_cent, cent_length] = find_centerline(P_start, P_end, DEM, R, ... 
                                                            'search_step',      loop_search_step, ... 
                                                            'search_angle',     loop_search_angle, ... 
                                                            'max_gradient',     loop_max_gradient, ... 
                                                            'window',           loop_window, ...
                                                            'max_length_factor',max_length_factor); 
        channel_length = sum(cent_length); % [pix]
        if ~isnan(channel_length)
            % found channel end, exit function
            disp(append("Found channel centerline on recursion no. ", string(it_no), ". "))
            return
        else
            % flip start and end points and try again
            it_no = it_no + 1; 
            [x_cent, y_cent, cent_length] = find_centerline(P_end, P_start, DEM, R, ... 
                                                                'search_step',      loop_search_step, ... 
                                                                'search_angle',     loop_search_angle, ... 
                                                                'max_gradient',     loop_max_gradient, ... 
                                                                'window',           loop_window, ...
                                                                'max_length_factor',max_length_factor); 
            channel_length = sum(cent_length); % [pix]
            if ~isnan(channel_length)
                % found channel end, flip and exit function
                x_cent = flip(x_cent); 
                y_cent = flip(y_cent); 
                cent_length = flip(cent_length); 
                disp(append("Found channel centerline on recursion no. ", string(it_no), " (reversed). "))
                return
            end
        end
    end
    
    disp("Warning: did not reach channel end (exceeded max. no. of centerline search recursions). ")
end

end % of function