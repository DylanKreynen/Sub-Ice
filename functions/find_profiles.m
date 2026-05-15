function [profiles, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, varargin)
%[profiles, x_prof, y_prof] = find_profiles(x_cent, y_cent, DEM, R, prof_length, prof_interval)
%Returns channel cross sectional profiles along channel centerline.
% basic idea:
% - find midpoint between two centerline locations (or position at user-defined interval)
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
% prof_length   = length of cross sectional profiles [m] (default: 2500m)
% prof_interval = spacing between profiles along centerline [m] (default: 0,
%                 meaning one profile per centerline segment at its midpoint)
%
% output:
% profiles = matrix containing profiles' sampled elevation [m]
% x_prof = matrix containing profiles' x coordinates [pix]
% Y_prof = matrix containing profiles' y coordinates [pix]
%
% (c) Dylan Kreynen
% University of Oslo
% 2024-2026

%% inputParser

% default parameter values
default_prof_length = 2500;
default_prof_interval = 0;

% parse input arguments
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validMapCellsRef = @(x) class(x) == "map.rasterref.MapCellsReference";
validxycent = @(x) (size(x, 1)==1 | size(x, 2)==1) && isnumeric(x(1)) && isscalar(x(1));
addRequired(p, 'x_cent', validxycent)
addRequired(p, 'y_cent', validxycent)
addRequired(p, 'DEM')
addRequired(p, 'R', validMapCellsRef)
addOptional(p, 'prof_length', default_prof_length, validScalarPosNum)
addParameter(p, 'prof_interval', default_prof_interval, validScalarPosNum)
parse(p, x_cent, y_cent, DEM, R, varargin{:});

prof_length   = p.Results.prof_length;
prof_interval = p.Results.prof_interval;


%% actual function

res = R.CellExtentInWorldX;        % resolution of DEM [m/pix]
DEM_int = griddedInterpolant(DEM); % interpolant for our DEM

% determine sampling step and number of pts on profile:
samp_step = res;
no_samp_pts = ceil(prof_length/samp_step);
no_samp_pts = 2*(ceil(no_samp_pts/2));

if prof_interval == 0
    % default: one profile per centerline segment, positioned at segment midpoint

    profiles = NaN(no_samp_pts+1, length(x_cent)-1);
    x_prof   = NaN(no_samp_pts+1, length(x_cent)-1);
    y_prof   = NaN(no_samp_pts+1, length(x_cent)-1);

    for i = 1:length(x_cent)-1
        P1 = [x_cent(i),   y_cent(i)  ];
        P2 = [x_cent(i+1), y_cent(i+1)];
        % construct sampling vector, perpendicular to and at section mid point
        [x_prof(:,i), y_prof(:,i)] = profile_query_pts(P1, P2, res, samp_step, no_samp_pts);
        % ATT! we might have to think about signs here, so profiles don't get flipped
        % profile observer should always be facing either up or downstream (seems OK, but which side is which?)
        profiles(:,i) = DEM_int(y_prof(:,i), x_prof(:,i));
    end

else
    % user-defined interval: distribute profiles at prof_interval spacing
    % along the centerline, independent of the centerline search step

    % cumulative length along centerline [pix]
    dx = diff(x_cent(:)');
    dy = diff(y_cent(:)');
    seg_lengths = sqrt(dx.^2 + dy.^2);
    cum_length  = [0, cumsum(seg_lengths)];

    % profile positions at prof_interval spacing [pix from centerline start]
    interval_pix = prof_interval / res;
    prof_s       = 0 : interval_pix : cum_length(end);
    no_profiles  = length(prof_s);

    profiles = NaN(no_samp_pts+1, no_profiles);
    x_prof   = NaN(no_samp_pts+1, no_profiles);
    y_prof   = NaN(no_samp_pts+1, no_profiles);

    for i = 1:no_profiles
        s = prof_s(i);

        % find which centerline segment this position falls in
        seg_idx = find(cum_length <= s, 1, 'last');
        seg_idx = min(seg_idx, length(x_cent)-1);  % clamp to last valid segment

        % interpolate position along segment
        frac  = (s - cum_length(seg_idx)) / seg_lengths(seg_idx);
        x_pos = x_cent(seg_idx) + frac * (x_cent(seg_idx+1) - x_cent(seg_idx));
        y_pos = y_cent(seg_idx) + frac * (y_cent(seg_idx+1) - y_cent(seg_idx));

        % build synthetic P1/P2 centered at profile position, aligned with the
        % local segment direction (profile_query_pts uses their midpoint and
        % the perpendicular of their connecting direction)
        unit_dx = (x_cent(seg_idx+1) - x_cent(seg_idx)) / seg_lengths(seg_idx);
        unit_dy = (y_cent(seg_idx+1) - y_cent(seg_idx)) / seg_lengths(seg_idx);
        P1 = [x_pos - 0.5*unit_dx, y_pos - 0.5*unit_dy];
        P2 = [x_pos + 0.5*unit_dx, y_pos + 0.5*unit_dy];

        [x_prof(:,i), y_prof(:,i)] = profile_query_pts(P1, P2, res, samp_step, no_samp_pts);
        profiles(:,i) = DEM_int(y_prof(:,i), x_prof(:,i));
    end

end

end

