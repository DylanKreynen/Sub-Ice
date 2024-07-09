function [x_samp, y_samp] = perp_prof_sampling(x_start, y_start, x_end, y_end, res, sampling_step, no_sampling_pts) 
%[x_samp, y_samp] = perp_prof_sampling(x_start, y_start, x_end, y_end, res, sampling_step, no_sampling_pts) 
%Returns a the coordinates of a sampling transect that is a certain
%distance away from a starting point (in the direction of an end point),
%perpendicular to that direction. 
%
% input: 
% x_start = x img coordinate start point [pix]
% y_start = y img coordinate start point [pix]
% x_end = x img coordinate end point [pix]
% y_end = y img coordinate end point [pix]
% res = image resolution [m/pixel]
% sampling_step = distance between sampling points on search profile [m]
% no_sampling_pts = number of sampling points on search profile [-]
%
% output: 
% x_samp = vector with x coordinates of sampling locations
% y_samp = vector with y coordinates of sampling locations
% 
% (c) Dylan Kreynen
% University of Oslo
% June - July 2024

sampling_step = sampling_step/res;  % from [m] to [pix]
no_sampling_pts = ceil(no_sampling_pts/2); % to end up with appropriate number of sampling profile points

x_mid = (x_start + x_end)/2; 
y_mid = (y_start + y_end)/2; 

sampling_constr = [linspace(-no_sampling_pts, 0, no_sampling_pts+1) linspace(1, no_sampling_pts, no_sampling_pts)]; 
sampling_constr = sampling_step*sampling_constr; 

% perpendicular direction: 
sample_dir = [y_end-y_start, x_start-x_end]; 
sample_dir = sample_dir./norm(sample_dir); 
% ATT! we might have to think about signs here, so profiles don't get flipped
% profile observer should always be facing either up or downstream (might already be OK)

% sampling x and y coords: 
x_samp = x_mid + sample_dir(1)*sampling_constr; 
y_samp = y_mid + sample_dir(2)*sampling_constr; 

end

