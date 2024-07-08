function [x_samp, y_samp] = profile_query_pts(P_start, P_end, res, samp_step, no_samp_pts) 
%Returns a the coordinates of a sampling transect that is a certain
%distance away from a starting point (in the direction of an end point),
%perpendicular to that direction. 
%
% input: 
% P_start = vector containing x and y coordinates of start point [pix coord]
% P_end = vector containing x and y coordinates of end point [pix coord]
% res = image resolution [m/pixel]
% samp_step = distance between sampling points on search profile [m]
% no_samp_pts = number of sampling points on search profile [-]
%
% output: 
% x_samp = vector with x coordinates of sampling locations
% y_samp = vector with y coordinates of sampling locations
% 
% (c) Dylan Kreynen
% University of Oslo
% June - July 2024

x_start = P_start(1); 
y_start = P_start(2); 
x_end = P_end(1); 
y_end = P_end(2); 

samp_step = samp_step/res;          % from [m] to [pix]
no_samp_pts = ceil(no_samp_pts/2);  % to end up with appropriate number of sampling profile points

x_mid = (x_start + x_end)/2; 
y_mid = (y_start + y_end)/2; 

sampling_constr = [linspace(-no_samp_pts, 0, no_samp_pts+1) linspace(1, no_samp_pts, no_samp_pts)]; 
sampling_constr = samp_step*sampling_constr; 

% perpendicular direction: 
sample_dir = [y_end-y_start, x_start-x_end]; 
sample_dir = sample_dir./norm(sample_dir); 
% ATT! we might have to think about signs here, so profiles don't get flipped
% profile observer should always be facing either up or downstream (might already be OK)

% sampling x and y coords: 
x_samp = x_mid + sample_dir(1)*sampling_constr; 
y_samp = y_mid + sample_dir(2)*sampling_constr; 

end

