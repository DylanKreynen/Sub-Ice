function [edge_idx, edge_coord, edge_elev] = find_edges(profiles, x_prof, y_prof, samp_step, slope_thr)
%Returns indices of channel cross sectional profiles that correspond with
%the channel's outer edges. Right now this is based on a slope threshold:
%returns the first indices where a slope threshold is met, áfter reaching 
%maximum slope. Right and left channel edges correspond to right and left 
%of channel centerline, looking downstream (from channel start to end). 
%
% input: 
% profiles = matrix containing profiles' sampled elevation [m]
% x_prof = matrix containing profiles' x coordinates [pix]
% y_prof = matrix containing profiles' y coordinates [pix]
% samp_step = distance between sampling points on profile [m]
% slope_thr = slope threshold for identifying channel edge [deg]
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
% June - July 2024

no_profs = size(profiles, 2);    % number of cross sectional profiles [-]
prof_length = size(profiles, 1); % length of cross sectional profiles [-] (no. of pts)

ledge_idx = NaN(no_profs, 1);   % left edge index ("upper" edge in previous versions)
redge_idx = NaN(no_profs, 1);   % right edge index ("lower" edge in previous versions)
lx = NaN(no_profs, 1);          % left edge x coord [pix]
ly = NaN(no_profs, 1);          % left edge y coord [pix]
lelev = NaN(no_profs, 1);       % channel elevation at left edge of profile [m]
rx = NaN(no_profs, 1);          % right edge x coord [pix]
ry = NaN(no_profs, 1);          % right edge y coord [pix]
relev = NaN(no_profs, 1);       % channel elevatoin at right edge of profile [m]
% TO DO: check whether this is still and always the case (left vs. right), 
% seems flipped now (for Venable). depends on start/end orientation
% left/right is correct when looking from END to START
% left/right is switched when looking from START to END

for i = 1:no_profs
    
    prof = profiles(:,i); 
    x_pr = x_prof(:,i); 
    y_pr = y_prof(:,i); 
    
    slope = gradient(prof, samp_step);  % slope in [rad]
    slope = abs(rad2deg(slope));        % now abs and in [deg]
    
    no_pts = ceil(prof_length/2);   % no of pts in "half" a profile (incl. midpoint) [-]
    
    % left channel edge (edge "to the left" of profile midpoint)
    % > find max slope to the left of channel midpoint 
    lslope = slope(1:no_pts); 
    [~, max_loc] = max(lslope);
    lslope = lslope(1:max_loc); 
    % > find first idx that satisfies threshold condition
    % (but after reaching maximum slope)
    idx = find(lslope < slope_thr); 
    if isempty(idx)       % if threshold is not reached
        idx = no_pts;     % take profile end as channel edge
    else
        idx = idx(end);   % only take the index closest to channel midpoint (after max slope)
    end 
    ledge_idx(i) = idx;
    
    % pixel coords and max elevation
    lx(i) = x_pr(idx); 
    ly(i) = y_pr(idx); 
    lelev(i) = prof(idx); 
    
    
    % right channel edge (edge "to the right" of profile midpoint)
    % > find max slope to the right of channel midpoint 
    rslope = flip(slope); % flipping profile slope for easier slicing
    rslope = rslope(1:no_pts); 
    [~, max_loc] = max(rslope);
    rslope = rslope(1:max_loc); 
    % > find first idx that satisfies threshold condition
    % (but after reaching maximum slope)
    idx = find(rslope < slope_thr); 
    if isempty(idx)       % if threshold is not reached
        idx = no_pts;     % take profile end as channel edge
    else
        idx = idx(end);   % only take the index closest to channel midpoint (after max slope)
    end 
    % profile slope was flipped! correcting for that:
    idx = prof_length+1 - idx; 
    redge_idx(i) = idx; 
    
    % pixel coords and max elevation
    rx(i) = x_pr(idx); 
    ry(i) = y_pr(idx); 
    relev(i) = prof(idx);
    
end

edge_idx = [ledge_idx redge_idx];  % [-]
edge_coord = [lx ly rx ry];        % [pix]
edge_elev = [lelev relev];         % [m]

