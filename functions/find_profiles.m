function [profiles, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, samp_step, no_samp_pts)
%Returns channel cross sectional profiles along channel centerline. 
% basic idea: 
% - find midpoint between two centerline locations 
% - find perpendicular direction to centerline
% - create sampling vector in this direction
% - sample DEM (interpolate) to find elevation profile
% - return profiles' coordinates and elevations
%
% input: 
% x_cent = vector with x coordinates of channel centerline [pix]
% y_cent = vector with y coordinates of channel centerline [pix]
% DEM = elevation data array [m] (use readgeoraster to read a geotiff)
% R = spatial referencing information for the array [-]
% samp_step = distance between sampling points on profile [m]
% no_samp_pts = number of sampling points on profile [-]
%
% output: 
% profiles = matrix containing profiles' sampled elevation [m]
% x_prof = matrix containing profiles' x coordinates [pix]
% Y_prof = matrix containing profiles' y coordinates [pix]
% 
% (c) Dylan Kreynen
% University of Oslo
% June - July 2024

res = R.CellExtentInWorldX; % resolution of DEM [m/pix]
x = 1:size(DEM, 1); 
y = 1:size(DEM, 2);
[X, Y] = meshgrid(y, x); 

% for storing profiles and coords: 
profiles = NaN(no_samp_pts+1, length(x_cent)-1);
x_prof = NaN(no_samp_pts+1, length(x_cent)-1);
y_prof = NaN(no_samp_pts+1, length(x_cent)-1);

for i = 1:length(x_cent)-1
    
    P1 = [x_cent(i), y_cent(i)];      % second most recent centerline point coords
    P2 = [x_cent(i+1), y_cent(i+1)];  % most recent centerline point coords
    
    % construct sampling vector, perpendicular to and at section mid point
    [x_prof(:,i), y_prof(:,i)] = profile_query_pts(P1, P2, res, samp_step, no_samp_pts);
    % ATT! we might have to think about signs here, so profiles don't get flipped
    % profile observer should always be facing either up or downstream (seems OK, but which side is which?)
    
    % sample (interpolate) elevation profile
    profiles(:,i) = interp2(X, Y, DEM, x_prof(:,i), y_prof(:,i));
    
end

end

