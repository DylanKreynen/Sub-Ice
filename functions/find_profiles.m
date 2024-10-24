function [profiles, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, varargin)
%[profiles, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, samp_step, no_samp_pts)
%Returns channel cross sectional profiles along channel centerline. 
% basic idea: 
% - find midpoint between two centerline locations 
% - find perpendicular direction to centerline
% - create sampling vector in this direction
% - sample DEM (interpolate) to find elevation profile
% - return profiles' coordinates and elevations
%
% required input: 
% x_cent = vector with x coordinates of channel centerline [pix]
% y_cent = vector with y coordinates of channel centerline [pix]
% DEM = elevation data array [m] (use readgeoraster to read a geotiff)
% R = spatial referencing information for the array [-]
%
% optional input: 
% samp_step = distance between sampling points on profile [m] (default: 100m)
% no_samp_pts = number of sampling points on profile [-] (default: 10)
%
% output: 
% profiles = matrix containing profiles' sampled elevation [m]
% x_prof = matrix containing profiles' x coordinates [pix]
% Y_prof = matrix containing profiles' y coordinates [pix]
% 
% (c) Dylan Kreynen
% University of Oslo
% June - July 2024

%% inputParser

% default parameter values
default_samp_step = 100; 
default_no_samp_pts =  10; 

% parse input arguments
p = inputParser; 
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validMapCellsRef = @(x) class(x) == "map.rasterref.MapCellsReference"; 
validxycent = @(x) (size(x, 1)==1 | size(x, 2)==1) && isnumeric(x(1)) && isscalar(x(1)); 
addRequired(p, 'x_cent', validxycent)
addRequired(p, 'y_cent', validxycent)
addRequired(p, 'DEM')
addRequired(p, 'R', validMapCellsRef)
addOptional(p, 'samp_step', default_samp_step, validScalarPosNum)
addOptional(p, 'no_samp_pts', default_no_samp_pts, validScalarPosNum)
parse(p, x_cent, y_cent, DEM, R, varargin{:}); 

samp_step = p.Results.samp_step; 
no_samp_pts = p.Results.no_samp_pts; 


%% actual function

res = R.CellExtentInWorldX;        % resolution of DEM [m/pix]
DEM_int = griddedInterpolant(DEM); % interpolant for our DEM

% for storing profiles and coords: 
no_samp_pts = 2*(ceil(no_samp_pts/2)); 
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
    %profiles(:,i) = interp2(DEM, x_prof(:,i), y_prof(:,i));
    profiles(:,i) = DEM_int(y_prof(:,i), x_prof(:,i)); 
    
end

end

