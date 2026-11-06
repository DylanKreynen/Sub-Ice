function T = channel_metrics(x_cent, y_cent, x_prof, y_prof, edge_idx, edge_coord, edge_elev, valid_edges, trough_elev, valid_trough, R, prof_interval)
%T = channel_metrics(x_cent, y_cent, x_prof, y_prof, edge_idx, edge_coord, edge_elev, valid_edges, trough_elev, valid_trough, R, prof_interval)
%Returns a table with channel metrics for each cross sectional profile of a
%channel: profile location, channel edge locations and elevations, trough
%elevation, channel width, half widths and depth.
%
% required input:
% x_cent, y_cent = vectors with centerline coordinates [pix] (see find_centerline())
% x_prof, y_prof = matrices with profiles' coordinates [pix] (see find_profiles())
% edge_idx, edge_coord, edge_elev, valid_edges = channel edges (see find_edges())
% trough_elev, valid_trough = channel trough (see find_trough())
% R = spatial referencing information for the DEM [-]
% prof_interval = spacing between profiles along centerline [m] (0 = one
%                 profile per centerline segment, at segment midpoint; see
%                 find_profiles())
%
% output:
% T = table, one row per profile, with columns:
%     prof_no           profile number [-]
%     dist_along_m      distance along centerline, from channel start [m]
%     x_center, y_center            profile center (centerline) [map units]
%     x_ledge, y_ledge, x_redge, y_redge   left and right edge [map units]
%     width_m           channel width (distance between edges) [m]
%     hw_left_m, hw_right_m   half width (distance of edge from centerline) [m]
%     elev_ledge_m, elev_redge_m   edge elevation [m]
%     elev_trough_m     trough elevation (profile center) [m]
%     depth_m           channel depth: elevation of left edge minus elevation
%                       of trough [m] (positive for a trough below the edge)
%     valid_ledge, valid_redge, valid_trough   outlier flags [-]
%                       (1 = original, 0 = outlier, left out or substituted)
%     (left and right as seen looking downstream, from start to end point)
%
% (c) Dylan Kreynen
% University of Oslo
% 2026


res = R.CellExtentInWorldX;         % resolution of DEM [m/pix]
no_profiles = size(x_prof, 2);
mid = (size(x_prof, 1) + 1)/2;      % index of profile center [-]

% distance along centerline of each profile (as placed by find_profiles())
seg_lengths = sqrt(diff(x_cent(:)).^2 + diff(y_cent(:)).^2)*res;   % [m]
if prof_interval == 0
    cum_length = [0; cumsum(seg_lengths)];
    dist_along = cum_length(1:end-1) + seg_lengths/2;
else
    dist_along = (0:no_profiles-1)'*prof_interval;
end

% profile centers and edges in map coordinates
[x_center, y_center] = intrinsicToWorld(R, x_prof(mid,:)', y_prof(mid,:)');
[x_ledge, y_ledge] = intrinsicToWorld(R, edge_coord(:,1), edge_coord(:,2));
[x_redge, y_redge] = intrinsicToWorld(R, edge_coord(:,3), edge_coord(:,4));

% half widths and width (left edge at higher profile indices than center)
hw_left = (edge_idx(:,1) - mid)*res;
hw_right = (mid - edge_idx(:,2))*res;
width = hw_left + hw_right;

% depth below left edge
depth = edge_elev(:,1) - trough_elev(:);

T = table((1:no_profiles)', dist_along, x_center, y_center, ...
    x_ledge, y_ledge, x_redge, y_redge, ...
    width, hw_left, hw_right, edge_elev(:,1), edge_elev(:,2), trough_elev(:), depth, ...
    valid_edges(:,1), valid_edges(:,2), valid_trough(:), ...
    'VariableNames', {'prof_no', 'dist_along_m', 'x_center', 'y_center', ...
    'x_ledge', 'y_ledge', 'x_redge', 'y_redge', ...
    'width_m', 'hw_left_m', 'hw_right_m', 'elev_ledge_m', 'elev_redge_m', 'elev_trough_m', 'depth_m', ...
    'valid_ledge', 'valid_redge', 'valid_trough'});

end
