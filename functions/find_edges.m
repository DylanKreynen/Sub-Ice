function [edge_idx, edge_coord, edge_elev, valid_edges, alongprof] = find_edges(profiles, x_prof, y_prof, res, varargin)
%[edge_idx, edge_coord, edge_elev] = find_edges(profiles, x_prof, y_prof,
%res, edge_method, slope_thr, slope_thr, min_width, max_width, peak_prom, sg_window, m_window)
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
% Near Peaks method is restricted to look within the minimum channel width
% and the maximum channel width. It detects if there are significant peaks
% within the area, and selects the closest peak to the channel. Otherwise
% it finds the knee point within the area.
%
% required input: 
% profiles = matrix containing profiles' sampled elevation [m]
% x_prof = matrix containing profiles' x coordinates [pix]
% y_prof = matrix containing profiles' y coordinates [pix]
% res = spatial resolution of DEM [m]
% 
% optional input: 
% edge_method = method to use to identify channel edges ("SlopeThreshold", "KneePoint" or "NearPeaks")
% knee_method = kneepoint algorithm to use ("LinearRegression" or "Kneedle", used with edge methods "KneePoint" and "NearPeaks") 
% slope_thr = slope threshold for identifying channel edge [deg] (default: 0 deg)
% min_width = minimum channel width [m] (default: 500m, set to 0 for no minimum width)
% max_width = maximum channel width [m] (set to 0 for maximum width = prof_length, default: 0)
% peak_prom = hMinPeakProm for findpeaks() [m], minimum prominence for channel edge (only used when edge_method = "NearPeaks")
% sg_window = window size for profile smoothing [m] (will be rounded, set to 0 for no smoothing)
% m_window = window size for edge smoothing [-] (no. of profile edges, set to 0 for no smoothing)
% keep_peaks = prevent peaks from being adjusted by along-channel edge smoothing (0 or 1, only used when edge_method = "NearPeaks")
% 
% optional input related to z-score outliers: 
% z_thr_elev =  z-score outlier threshold, elevation of edge (set to 0 to skip outlier identification, default: 0)
%               optional: [left_thr right_thr] to use different thresholds for the left and right edge
% z_thr_idx =   z-score outlier threshold, profile index of edge (set to 0 to skip outlier identification, default: 0)
%               optional: [left_thr right_thr] to use different thresholds for the left and right edge
% edge_subst_window = window size for outlier substitution (moving median filter, set to 0 to leave out outliers altogether, default: 5)  
%
% output: 
% edge_idx = matrix containing profile indices corr. to channel edges [-]
%               1st col: idx w.r.t. left channel edge
%               2nd col: idx w.r.t. right channel edge
% edge_coord = matrix containing channel edge coordinates [pix]
%               1st col: left channel edge x coordinate
%               2nd col: left channel edge y coordinate
%               3rd col: right channel edge x coordinate
%               4th col: right channel edge y coordinate
% edge_elev = matrix containing channel edge elevations [m]
%               1st col: left channel edge elevation
%               2nd col: right channel edge elevation]
% valid_edges = matrix containing an edge outlier flag [-]
%               1 = valid/original channel edge; 0 = outlier (left out or substituted)
%               1st col: left channel edge outlier flag
%               2nd col: right channel edge outlier flag
% 
% (c) Dylan Kreynen
% University of Oslo
% 2024-2026
% 
% "NearPeaks" method by Adam Isaak
% University of Manitoba
% 2026



%% input parser

% default parameter values
default_slope_thr = 0; 
default_min_width = 500; 
default_max_width = 0;
default_peak_prom = 1;
default_sg_window = 0; 
default_m_window = 0; 
default_edge_method = "KneePoint";
default_keep_pks = false;
default_knee_method = "LinearRegression"; 
default_z_thr_elev = 0; 
default_z_thr_idx = 0; 
default_edge_subst_window = 5;

% parse input arguments
p = inputParser; 
validScalarPosNum = @(x) isnumeric(x) && isscalar(x); % && (x >= 0);
validMaxMinWidths = @(x) (isvector(x) && all(x(:) >= 0) && length(x) == 2) || (isnumeric(x) && isscalar(x) && (x >= 0));
validEdgeMethod = @(x) convertCharsToStrings(x)=="SlopeThreshold" | convertCharsToStrings(x)=="KneePoint" | convertCharsToStrings(x)=="NearPeaks";
validKneeMethod = @(x) convertCharsToStrings(x)=="LinearRegression" | convertCharsToStrings(x)=="Kneedle";
addRequired(p, 'profiles')
addRequired(p, 'x_prof')
addRequired(p, 'y_prof')
addRequired(p, 'res', validScalarPosNum)
addOptional(p, 'slope_thr', default_slope_thr, validScalarPosNum)
addOptional(p, 'min_width', default_min_width, validMaxMinWidths)
addOptional(p, 'max_width', default_max_width, validMaxMinWidths)
addOptional(p, 'peak_prom', default_peak_prom, validScalarPosNum)
addOptional(p, 'sg_window', default_sg_window, validScalarPosNum)
addOptional(p, 'm_window', default_m_window, validScalarPosNum)
addOptional(p, 'edge_method', default_edge_method, validEdgeMethod)
addOptional(p, 'keep_pks', default_keep_pks, validScalarPosNum)
addOptional(p, 'knee_method', default_knee_method, validKneeMethod)
addOptional(p, 'z_thr_elev', default_z_thr_elev, validMaxMinWidths)
addOptional(p, 'z_thr_idx', default_z_thr_idx, validMaxMinWidths)
addOptional(p, 'edge_subst_window', default_edge_subst_window, validScalarPosNum)
parse(p, profiles, x_prof, y_prof, res, varargin{:}); 

slope_thr = p.Results.slope_thr; 
min_width = p.Results.min_width; 
max_width = p.Results.max_width;
peak_prom = p.Results.peak_prom;
sg_window = p.Results.sg_window; 
m_window = p.Results.m_window; 
edge_method = convertCharsToStrings(p.Results.edge_method);
keep_pks = p.Results.keep_pks;
knee_method = p.Results.knee_method;
z_thr_elev = p.Results.z_thr_elev; 
z_thr_idx = p.Results.z_thr_idx; 
edge_subst_window = p.Results.edge_subst_window; 


%% actual function

no_profs = size(profiles, 2);    % number of cross sectional profiles [-]
prof_length = size(profiles, 1); % length of cross sectional profiles [-] (no. of pts)
samp_step = sqrt((x_prof(1,1)-x_prof(2,1))^2 + (y_prof(1,1)-y_prof(2,1))^2); % [pix]
samp_step = samp_step*res;       % distance between sampling points, now in [m]

% channel edge should be at least this distance away from channel centerline


if ~isscalar(max_width)
    lmax_width = ceil(max_width(1)/samp_step);
    rmax_width = ceil(max_width(2)/samp_step);
else
    max_width = max_width/samp_step;    % from [m] to [-] (index)
    max_width = ceil(max_width/2);      % half-distance
    lmax_width = max_width;
    rmax_width = max_width;
end

% max_width = 0 means no limit: use full half-profile length on each side
if lmax_width == 0; lmax_width = floor(prof_length/2); end
if rmax_width == 0; rmax_width = floor(prof_length/2); end

if ~isscalar(min_width)
    lmin_width = ceil(min_width(1)/samp_step);
    rmin_width = ceil(min_width(2)/samp_step);
else
    min_width = min_width/samp_step;    % from [m] to [-] (index)
    min_width = ceil(min_width/2);      % half-distance
    lmin_width = min_width;
    rmin_width = min_width;
end

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
relev = NaN(no_profs, 1);        % channel elevation at right edge of profile [m]
% left/right is correct when looking from START to END

% smoothing array
ledge_sm = true(no_profs, 1);
redge_sm = true(no_profs, 1);

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
    % note: should we only smooth for slope threshold method, or also kneepoint/nearpeaks?

    if edge_method == "SlopeThreshold"
        % derivative to find slope
        slope = gradient(prof, samp_step);  % slope in [rad]
        
        % right channel edge (edge "to the right" of profile midpoint)
        rslope = slope(1:no_pts-min_width); 
    
        % find first idx that satisfies threshold condition
        idx = find(rslope > -slope_thr); 
        if isempty(idx)                     % if threshold is not reached
            idx = 1;                        % take profile end as channel edge
        else
            idx = idx(end);                 % only take the index closest to channel midpoint (after max slope)
        end 
        redge_idx(i) = idx;
        % clamp to min/max width
        if rmin_width > 0; redge_idx(i) = min(redge_idx(i), no_pts - rmin_width); end
        redge_idx(i) = max(redge_idx(i), no_pts - rmax_width);

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
        idx = prof_length - idx;
        ledge_idx(i) = idx;
        % clamp to min/max width
        if lmin_width > 0; ledge_idx(i) = max(ledge_idx(i), no_pts + lmin_width); end
        ledge_idx(i) = min(ledge_idx(i), no_pts + lmax_width);

    elseif edge_method == "KneePoint"

        % right channel edge
        rprof = prof(1:no_pts);
        rprof = flip(rprof);
        [~, idx] = knee_pt(rprof, 'knee_method', knee_method);
        % profile was flipped, correcting for that:
        redge_idx(i) = no_pts - idx;
        % clamp to min/max width
        if rmin_width > 0; redge_idx(i) = min(redge_idx(i), no_pts - rmin_width); end
        redge_idx(i) = max(redge_idx(i), no_pts - rmax_width);

        % left channel edge
        lprof = prof(no_pts:end);
        [~, idx] = knee_pt(lprof, 'knee_method', knee_method);
        ledge_idx(i) = no_pts + idx;
        % clamp to min/max width
        if lmin_width > 0; ledge_idx(i) = max(ledge_idx(i), no_pts + lmin_width); end
        ledge_idx(i) = min(ledge_idx(i), no_pts + lmax_width);

    elseif edge_method == "NearPeaks"
        % right channel edge
        out_th = no_pts-rmax_width;      % index of outer threshold

        rprof = prof(1:no_pts-rmin_width);

        % find the peaks along the right channel edge
        [~, pk] = findpeaks(rprof((no_pts-rmax_width):end), MinPeakProminence=peak_prom);

        if isempty(pk)                                  % if no peaks are found
            [~, idx] = knee_pt(rprof(out_th:end), 'knee_method', knee_method);      % find the knee point in the search area
            redge_idx(i) = idx+out_th;                  
            redge_sm(i) = true;                         % set as filterable
        else
            idx = ceil(pk(end)+out_th);                 % find the index of the peak
            redge_idx(i) = idx;                         
            redge_sm(i) = false;                        % set to preserve during filtering
        end 

        % left channel edge
        out_th = no_pts-lmax_width;      % index of outer threshold
        
        lprof = flip(prof); 
        lprof = lprof(1:no_pts-lmin_width);

        % find the peaks along the left channel edge
        [~, pk] = findpeaks(lprof(out_th:end), MinPeakProminence=peak_prom);

        if isempty(pk)                                  % if no peaks are found
            [~, idx] = knee_pt(lprof(out_th:end), 'knee_method', knee_method);      % find the knee point in the search area
            idx = prof_length - (idx+out_th);         % profile was flipped! correcting for that:
            ledge_sm(i) = true;                         % set as filterable
        else                                            
            idx = ceil(pk(end)+out_th);                 % find the index of the peak
            idx = prof_length - idx;                  % profile was flipped! correcting for that:
            ledge_sm(i) = false;                        % set to preserve during filtering
        end 

        ledge_idx(i) = idx;

    else
        error("Invalid edge method. Check find_edges() parameters, set edge_method to 'SlopeThreshold', 'KneePoint' or 'NearPeaks'.")
    end
end


% smooth edges using median filter (if window ~= 0) 
if m_window ~= 0
    ledge_idx_filt = ceil(medfilt1(ledge_idx, m_window, [], 1, 'truncate')); 
    redge_idx_filt = ceil(medfilt1(redge_idx, m_window, [], 1, 'truncate'));


    if keep_pks == 1         % preserve the peaks position from smoothing by resetting them back
        for i =1:no_profs
            if ~ledge_sm(i)
                ledge_idx_filt(i) = ledge_idx(i);
            end
    
            if ~redge_sm(i)
                redge_idx_filt(i) = redge_idx(i);
            end
        end
    end
else
    ledge_idx_filt = ledge_idx; 
    redge_idx_filt = redge_idx; 
end

% note: now we smooth the edges using a median filter, but take the exact
% elevation at the (smoothed) edge coordinates - probably not optimal
% (z-score outlier filtering, down below, probably better implemented)


%% edge filtering and output preparation
%  identify and deal with outliers, based on elevation and width
%  for now: z-score outliers and moving median substitution

% deal with z-score thresholds, might be different for left and right edge
if length(z_thr_elev) == 2
    z_thr_lelev = z_thr_elev(1); 
    z_thr_relev = z_thr_elev(2); 
elseif length(z_thr_elev) == 1
    z_thr_lelev = z_thr_elev; 
    z_thr_relev = z_thr_elev;
else
    error("Check dimensions of z_thr_elev. Should be 1x1 or 1x2. ")
end

if length(z_thr_idx) == 2
    z_thr_lidx = z_thr_idx(1);
    z_thr_ridx = z_thr_idx(2); 
elseif length(z_thr_idx) == 1; 
    z_thr_lidx = z_thr_idx;
    z_thr_ridx = z_thr_idx; 
else
    error("Check dimensions of z_thr_idx. Should be 1x1 or 1x2. ")
end

% raise warning when median smoothing ánd outlier filtering
if any([z_thr_lelev, z_thr_relev, z_thr_lidx, z_thr_ridx]) && m_window ~= 0
    warning("Note: performing along-channel edge smoothing before outlier filtering. Suggestion: set m_window = 0 to skip along-channel smoothing of channel edges.")
end

% fetch edge elevations
for i = 1:no_profs
    prof = profiles(:,i); 
    lelev(i) = prof(ledge_idx_filt(i));
    relev(i) = prof(redge_idx_filt(i)); 
end

% z-score for edge profile index
mean_idx = mean(ledge_idx_filt); 
std_idx = std(ledge_idx_filt); 
z_lidx = (ledge_idx_filt-mean_idx)./std_idx; % left edge
mean_idx = mean(redge_idx_filt); 
std_idx = std(redge_idx_filt); 
z_ridx = (redge_idx_filt-mean_idx)./std_idx; % right edge

% z-score for edge elevation
mean_elev = mean(lelev); 
std_elev = std(lelev); 
z_lelev = (lelev-mean_elev)./std_elev; % left edge
mean_elev = mean(relev); 
std_elev = std(relev); 
z_relev = (relev-mean_elev)./std_elev; % right edge

% flag outliers (1 = keep, 0 = outlier)
% threshold = 0 means "skip filtering" for that criterion
valid_ledge = ones(size(ledge_idx_filt));
if z_thr_lidx > 0;  valid_ledge(abs(z_lidx)  > z_thr_lidx)  = 0; end
if z_thr_lelev > 0; valid_ledge(abs(z_lelev) > z_thr_lelev) = 0; end
valid_redge = ones(size(redge_idx_filt));
if z_thr_ridx > 0;  valid_redge(abs(z_ridx)  > z_thr_ridx)  = 0; end
if z_thr_relev > 0; valid_redge(abs(z_relev) > z_thr_relev) = 0; end
valid_edges = [valid_ledge valid_redge]; % both edges

% substitute outliers
% two options: 
%   - leave out outliers (replace with NaN)
%   - substitute values based on moving median filter
% temporarily set filtered out edges to NaN

lelev(valid_ledge==0) = NaN; 
relev(valid_redge==0) = NaN; 
ledge_idx_filt(valid_ledge==0) = NaN; 
redge_idx_filt(valid_redge==0) = NaN; 

if edge_subst_window ~= 0 
% use moving median filter to replace outliers

    % construct moving median substitute values
    med_lelev = movmedian(lelev, edge_subst_window, "omitnan"); 
    med_relev = movmedian(relev, edge_subst_window, "omitnan"); 
    med_lidx = movmedian(ledge_idx_filt, edge_subst_window, "omitnan"); 
    med_ridx = movmedian(redge_idx_filt, edge_subst_window, "omitnan"); 
    
    % substitute where necessary
    lelev(valid_ledge==0) = med_lelev(valid_ledge==0); 
    relev(valid_redge==0) = med_relev(valid_redge==0); 
    ledge_idx_filt(valid_ledge==0) = med_lidx(valid_ledge==0); 
    redge_idx_filt(valid_redge==0) = med_ridx(valid_redge==0); 

end

% convert to integers for indexing, prep output
% use round() instead of int32() to preserve NaN for filtered-out outliers
ledge_idx_filt = round(ledge_idx_filt);
redge_idx_filt = round(redge_idx_filt);
edge_idx = [ledge_idx_filt redge_idx_filt];
edge_elev = [lelev relev]; 

% fetch edge image coordinates

for i = 1:no_profs

    x_pr = x_prof(:,i); 
    y_pr = y_prof(:,i); 

    lidx = ledge_idx_filt(i);
    ridx = redge_idx_filt(i);

    % image coordinates [pix]

    if isnan(lidx) || ~lidx
        lx(i) = NaN;
        ly(i) = NaN;
    else
        lx(i) = x_pr(lidx); 
        ly(i) = y_pr(lidx);
    end

    if isnan(ridx) || ~ridx
        rx(i) = NaN; 
        ry(i) = NaN; 
    else
        rx(i) = x_pr(ridx); 
        ry(i) = y_pr(ridx); 
    end

end

edge_coord = [lx ly rx ry];        % [pix]