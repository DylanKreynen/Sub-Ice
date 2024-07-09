function [x_samp, y_samp] = perp_search_sampling(x_start, y_start, x_end, y_end, res, search_step, sampling_step, no_sampling_pts, first_after_start)
%[x_samp, y_samp] = perp_search_sampling(x_start, y_start, x_end, y_end, res, search_step, sampling_step, no_sampling_pts, first_after_start)
%Returns the coordinates of a sampling transect that is a certain
%distance away from a starting point (in the direction of an end point),
%perpendicular to that direction. 
%
% input: 
% x_start = x img coordinate start point [pix]
% y_start = y img coordinate start point [pix]
% x_end = x img coordinate end point [pix]
% y_end = y img coordinate end point [pix]
% res = image resolution [m/pixel]
% search_step = distance to step away from P_start to construct search profile [m]
% sampling_step = distance between sampling points on search profile [m]
% no_sampling_pts = number of sampling points on search profile [-]
% first_after_start = 1 or 0, whether we're looking for the first point 
% after channel start point
%
% output: 
% x_samp = vector with x coordinates of sampling locations
% y_samp = vector with y coordinates of sampling locations
% 
% (c) Dylan Kreynen
% University of Oslo
% June - July 2024

search_step = search_step/res;      % from [m] to [pix]
sampling_step = sampling_step/res;  % from [m] to [pix]

no_sampling_pts = ceil(no_sampling_pts/2);      % to end up with appropriate number of sampling profile points

search_dir = [x_end-x_start, y_end-y_start];    % perpendicular direction to start-end vector (search direction)
search_dir = search_dir./norm(search_dir);      % now a unit vector

% midpoint one search step away: 
if first_after_start == 1
    sampling_mid_x = x_start + search_step*search_dir(1); 
    sampling_mid_y = y_start + search_step*search_dir(2); 
else
    sampling_mid_x = x_end + search_step*search_dir(1); 
    sampling_mid_y = y_end + search_step*search_dir(2); 
end

sampling_constr = [linspace(-no_sampling_pts, 0, no_sampling_pts+1) linspace(1, no_sampling_pts, no_sampling_pts)]; 
sampling_constr = sampling_step*sampling_constr; 

% perpendicular direction: 
sample_dir = [y_end-y_start, x_start-x_end]; 
sample_dir = sample_dir./norm(sample_dir); 

% x and y coords of sampling points: 
x_samp = sampling_mid_x + sample_dir(1)*sampling_constr; 
y_samp = sampling_mid_y + sample_dir(2)*sampling_constr; 

% idea: implement the same idea for a semi-circle sampling vector, rather
% than a straight-line section - that way all sampling points are an equal
% distance away from the last known channel centerline point

% % attempt at semicircle sampling vector: 
% search_angle = 90; % [deg]
% search_angle = deg2rad(search_angle); % [rad]
% search_start = atan(search_dir(2)/search_dir(1)) + search_angle/2;   % in [rad]
% search_stop = atan(search_dir(2)/search_dir(1)) - search_angle/2;    % in [rad]
% th = linspace(search_start, search_stop, no_search_pts);
% x = search_step*cos(th) + x_start;
% y = search_step*sin(th) + y_start;
% % (kinda works but not in quite the right direction)

end

