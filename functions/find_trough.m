function [trough_elev, trough_depth, valid_trough] = find_trough(profiles, edge_elev, varargin)
%[trough_elev, trough_depth, valid_trough] = find_trough(profiles, edge_elev)
%
% Returns the absolute elevation and depth (w.r.t. left channel edge) of
% the channel trough. Channel trough is assumed to be located at exactly
% the center of the cross sectional profiles. Outliers can be flagged
% (based on a threshold on z-score) and left out or replaced using a moving
% mean filter.

% required input:
% profiles = matrix containing profiles' sampled elevation [m]
% edge_elev = matrix containing channel edge elevations [m]
%               1st col: left channel edge elevation
%               2nd col: right channel edge elevation
%
% optional input:
% z_thr = z-score outlier threshold, elevation of trough (set to 0 to skip outlier identification, default: 0)
% subs_window = window for outlier substitution (moving mean filter, set to 0 to leave out outliers altogether, default: 3)
%
% output:
% trough_elev = vector containing absolute trough elevation [m]
% trough_depth = vector containing trough depth w.r.t. left channel edge [m]
% valid_trough = vector containing a trough outlier flag [-]
%               1 = valid/original trough; 0 = outlier (left out or substituted)
%
% (c) Dylan Kreynen
% University of Oslo
% 2024-2026


%% input parser

default_z_thr = 0;
default_subst_window = 3;

p = inputParser;
validScalarNonNeg = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addRequired(p, 'profiles')
addRequired(p, 'edge_elev')
addOptional(p, 'z_thr', default_z_thr, validScalarNonNeg)
addOptional(p, 'subst_window', default_subst_window, validScalarNonNeg)
parse(p, profiles, edge_elev, varargin{:});

z_thr = p.Results.z_thr;
subst_window = p.Results.subst_window;


%% actual function

no_profs = size(profiles, 2);
prof_length = size(profiles, 1);
no_pts = ceil(prof_length / 2);  % center index = trough location

% extract trough elevation at profile center
trough_elev = profiles(no_pts, :)';  % [no_profs x 1]

% z-score outlier detection on trough elevation
valid_trough = ones(no_profs, 1);

if z_thr > 0
    mean_elev = mean(trough_elev, 'omitnan');
    std_elev = std(trough_elev, 'omitnan');
    z_scores = (trough_elev - mean_elev) ./ std_elev;
    valid_trough(abs(z_scores) > z_thr) = 0;
end

% set outliers to NaN
trough_elev(valid_trough == 0) = NaN;

% substitute outliers using moving mean filter
if subst_window ~= 0
    mean_trough = movmean(trough_elev, subst_window, 'omitnan');
    trough_elev(valid_trough == 0) = mean_trough(valid_trough == 0);
end

% compute depth w.r.t. left channel edge (negative when trough is below edge)
trough_depth = trough_elev - edge_elev(:, 1);

end
