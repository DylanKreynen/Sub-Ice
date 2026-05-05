function [edge_idx, edge_coord, edge_elev, alongprof] = find_edges(profiles, x_prof, y_prof, res, varargin)
%[edge_idx, edge_coord, edge_elev] = find_edges(profiles, x_prof, y_prof,
%res, edge_method, slope_thr, slope_thr, min_width, sg_window, m_window)
%
% Returns indices of channel cross sectional profiles that correspond with
% the channel's outer edges, based on either a slope threshold or knee point
% method. Edges can be a certain min. distance away from channel center
% line. Right and left channel edges correspond to right and left of channel 
% centerline, looking downstream (from channel start to end). 
% 
% Profiles can be smoothed across-channel using a Savitsky-Golay filter 
% before applying edge detection (especially useful for the slope threshold
% method). Edge coordinates can be smoothed along-profile using an optional 
% median filter (applicable to all methods). This version of this function 
% also outputs "along profiles" for output figures (work in progress).  
%
% required input: 
% profiles = matrix containing profiles' sampled elevation [m]
% x_prof = matrix containing profiles' x coordinates [pix]
% y_prof = matrix containing profiles' y coordinates [pix]
% res = spatial resolution of DEM [m]
% 
% optional input: 
% edge_method = method to use to find channel edges (default "SlopeThreshold" alternatively "KneePoint")
% slope_thr = slope threshold for identifying channel edge [deg] (default: 0 deg)
% min_width = minimum channel width (edge must be half min. width away from center line) [m] (default: 500m)
% sg_window = window size for profile smoothing [m] (will be rounded, set to 0 for no smoothing)
% m_window = window size for edge smoothing [-] (no. of profile edges, set to 0 for no smoothing)
%
% output: 
% edge_idx = matrix containing profile indices corr. to channel edges [-]
%               1st col: idx w.r.t. left channel edge
%               2nd col: idx w.r.t. right channel edge
% egde_coord = matrix containing channel edge coordinates [pix]
%               1st col: left channel edge x coordinate
%               2nd col: left channel edge y coordinate
%               3rd col: right channel edge x coordinate
%               4th col: right channel edge y coordinate
% edge_elev = matrix containing channel edge elevations [m]
%               1st col: left channel edge elevation
%               2nd col: right channel edge elevation]
% 
% (c) Dylan Kreynen
% University of Oslo
% 2024-2026



%% input parser

% default parameter values
default_slope_thr = 0; 
default_min_width = 500; 
default_sg_window = 0; 
default_m_window = 0; 
default_edge_method = "KneePoint"; 

% parse input arguments
p = inputParser; 
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validEdgeMethod = @(x) convertCharsToStrings(x)=="SlopeThreshold" | convertCharsToStrings(x)=="KneePoint"; 
addRequired(p, 'profiles')
addRequired(p, 'x_prof')
addRequired(p, 'y_prof')
addRequired(p, 'res', validScalarPosNum)
addOptional(p, 'slope_thr', default_slope_thr, validScalarPosNum)
addOptional(p, 'min_width', default_min_width, validScalarPosNum)
addOptional(p, 'sg_window', default_sg_window, validScalarPosNum)
addOptional(p, 'm_window', default_m_window, validScalarPosNum)
addOptional(p, 'edge_method', default_edge_method, validEdgeMethod)
parse(p, profiles, x_prof, y_prof, res, varargin{:}); 

slope_thr = p.Results.slope_thr; 
min_width = p.Results.min_width; 
sg_window = p.Results.sg_window; 
m_window = p.Results.m_window; 
edge_method = convertCharsToStrings(p.Results.edge_method); 


%% actual function

no_profs = size(profiles, 2);    % number of cross sectional profiles [-]
prof_length = size(profiles, 1); % length of cross sectional profiles [-] (no. of pts)
samp_step = sqrt((x_prof(1,1)-x_prof(2,1))^2 + (y_prof(1,1)-y_prof(2,1))^2); % [pix]
samp_step = samp_step*res;       % distance between sampling points, now in [m]

% channel edge should be at least this distance away from channel centerline
min_width = min_width/samp_step;    % from [m] to [-] (index)
min_width = ceil(min_width/2);      % half-distance

slope_thr = deg2rad(slope_thr);         % from [deg] to [rad]
sg_window = ceil(sg_window/samp_step);  % from [m] yo [-] (index)
if sg_window ~= 0 && mod(sg_window,2) == 0
    % make window odd
    sg_window = sg_window + 1; 
end

ledge_idx = NaN(no_profs, 1);    % left edge index ("upper" edge in previous versions)
redge_idx = NaN(no_profs, 1);    % right edge index ("lower" edge in previous versions)
lx = NaN(no_profs, 1);           % left edge x coord [pix]
ly = NaN(no_profs, 1);           % left edge y coord [pix]
lelev = NaN(no_profs, 1);        % channel elevation at left edge of profile [m]
rx = NaN(no_profs, 1);           % right edge x coord [pix]
ry = NaN(no_profs, 1);           % right edge y coord [pix]
relev = NaN(no_profs, 1);        % channel elevatoin at right edge of profile [m]
% left/right is correct when looking from START to END

% along profiles
% TO DO: now hard coded, better would be based on channel width (?)
alongprofidx = [250 500 750 1000 1250 1500]; % [m]
alongprofidx = floor(alongprofidx./samp_step); 
alongprof = NaN(no_profs, 2*length(alongprofidx)+1); 

for i = 1:no_profs

    prof = profiles(:,i); 

    no_pts = ceil(prof_length/2);                       % no of pts in "half" a profile (incl. midpoint) [-]

    % along profiles (no smoothing)
    ralongprof = prof([no_pts no_pts+alongprofidx]);    % incl. "centerline"
    lalongprof = prof(fliplr(alongprofidx));            % excl. "centerline"
    alongprof(i,:) = [lalongprof' ralongprof']; 

    % smooth profile using Savitzky-Golay filter 
    if sg_window ~= 0
        prof = sgolayfilt(prof, 3, sg_window); 
    end
    % note: should we only smooth for slope threshold method, or also knee?

    if edge_method == "SlopeThreshold"
        % derivative to find slope
        slope = gradient(prof, samp_step);  % slope in [rad]
        
        % right channel edge (edge "to the right" of profile midpoint)
        rslope = slope(1:no_pts-min_width); 
    
        % find first idx that satisfies threshold condition
        idx = find(rslope > slope_thr); 
        if isempty(idx)                     % if threshold is not reached
            idx = 1;                        % take profile end as channel edge
        else
            idx = idx(end);                 % only take the index closest to channel midpoint (after max slope)
        end 
        redge_idx(i) = idx;
    
        % left channel edge (edge "to the left" of profile midpoint)
        lslope = flip(slope);               % flipping profile slope for easier slicing
        lslope = lslope(1:no_pts-min_width); 
    
        % find first idx that satisfies threshold condition
        idx = find(lslope < slope_thr); 
        if isempty(idx)                     % if threshold is not reached
            idx = 1;                        % take profile end as channel edge
        else
            idx = idx(end);                 % only take the index closest to channel midpoint (after max slope)
        end 

        % profile slope was flipped! correcting for that:
        idx = prof_length+1 - idx; 
        ledge_idx(i) = idx; 

    elseif edge_method == "KneePoint"
        % right channel edge
        rprof = prof(1:no_pts-min_width); 
        % rprof = prof(1:no_pts); 
        [~, idx] = knee_pt(rprof); 
        redge_idx(i) = idx;

        % left channel edge
        lprof = flip(prof); 
        lprof = lprof(1:no_pts-min_width);
        % lprof = lprof(1:no_pts);
        [~, idx] = knee_pt(lprof); 

        % profile was flipped! correcting for that:
        idx = prof_length+1 - idx; 
        ledge_idx(i) = idx; 

    else
        error("Invalid edge method. Check find_edges() parameters, set edge_method to 'SlopeTreshold' or 'KneePoint'.")
    end
end


% smooth edges using median filter (if window ~= 0) 
if m_window ~= 0
    ledge_idx_filt = ceil(medfilt1(ledge_idx, m_window, [], 1, 'truncate')); 
    redge_idx_filt = ceil(medfilt1(redge_idx, m_window, [], 1, 'truncate')); 
else
    ledge_idx_filt = ledge_idx; 
    redge_idx_filt = redge_idx; 
end
edge_idx = [ledge_idx_filt redge_idx_filt];  % [-]


% go from indices to coordinates and fetch edge elevations

for i = 1:no_profs

    prof = profiles(:,i); 
    x_pr = x_prof(:,i); 
    y_pr = y_prof(:,i); 

    lidx = ledge_idx_filt(i);
    ridx = redge_idx_filt(i);

    % pixel coords and elevation
    rx(i) = x_pr(ridx); 
    ry(i) = y_pr(ridx); 
    relev(i) = prof(ridx); 

    % pixel coords and elevation
    lx(i) = x_pr(lidx); 
    ly(i) = y_pr(lidx); 
    lelev(i) = prof(lidx);

end

edge_coord = [lx ly rx ry];        % [pix]
edge_elev = [lelev relev];         % [m]

% note: now we smooth the edges using a median filter, but take the exact
% elevation at the (smoothed) edge coordinates - probably not optimal