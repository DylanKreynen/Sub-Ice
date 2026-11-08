function plot_channel_atlas_page(f, DEM, R, x_cent, y_cent, edge_coord, x_prof, y_prof, profiles, edge_idx, edge_elev, T, label, varargin)
%plot_channel_atlas_page(f, DEM, R, x_cent, y_cent, edge_coord, x_prof, y_prof, profiles, edge_idx, edge_elev, T, label, ...)
%Plots a channel atlas page: an overview of one channel on a single figure,
%in the shared Sub-Ice figure style (see fig_style()), with four panels:
%   (a) map of the channel (see plot_overview_map())
%   (b) cross-sectional profiles (see plot_profile_figure())
%   (c) elevation of the trough and channel edges along the channel
%   (d) channel width and depth below the left edge along the channel
%
% required input:
% f = figure handle (the figure is cleared)
% DEM = elevation data array [m] (no data values set to NaN)
% R = spatial referencing information for the DEM [-]
% x_cent, y_cent = centerline coordinates [pix] (see find_centerline())
% edge_coord = channel edge coordinates [pix] (see find_edges())
% x_prof, y_prof = profile coordinates [pix] (see find_profiles())
% profiles = profiles' sampled elevation [m] (see find_profiles())
% edge_idx, edge_elev = channel edges (see find_edges())
% T = table with channel metrics per profile (see channel_metrics())
% label = channel label, used as page title
%
% optional input:
% subtitle = page subtitle, e.g. the edge detection method (default: none)
% clims, hillshade, transects, bridge_gaps = map options, see plot_overview_map()
% min_width, max_width = edge search limits [m], see plot_profile_figure()
%
% (c) Dylan Kreynen
% University of Oslo
% 2026


%% inputParser

p = inputParser;
addParameter(p, 'subtitle', "")
addParameter(p, 'clims', [])
addParameter(p, 'hillshade', true)
addParameter(p, 'transects', false)
addParameter(p, 'bridge_gaps', true)
addParameter(p, 'min_width', 0)
addParameter(p, 'max_width', 0)
parse(p, varargin{:});
o = p.Results;


%% actual function

S = fig_style();
res = R.CellExtentInWorldX;
clf(f)
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, string(label), 'FontWeight', 'bold', 'FontSize', S.font_title + 1, 'Color', S.ink, 'Interpreter', 'none')
if strlength(string(o.subtitle)) > 0
    subtitle(tl, o.subtitle, 'FontSize', S.font_small + 0.5, 'Color', S.muted, 'Interpreter', 'none')
end

% (a) map of this channel only
ax = nexttile(tl, 1);
plot_overview_map(ax, DEM, R, {x_cent}, {y_cent}, {edge_coord}, {x_prof}, {y_prof}, string(label), ...
    'clims', o.clims, 'hillshade', o.hillshade, 'transects', o.transects, 'bridge_gaps', o.bridge_gaps, ...
    'title', "(a)  Map")

% (b) cross-sectional profiles, absolute elevation
ax = nexttile(tl, 2);
prof_dist = ((1:size(profiles, 1)) - (size(profiles, 1)+1)/2)*res;
plot_profile_figure(ax, prof_dist, profiles, edge_idx, edge_elev, ...
    'title', "(b)  Cross-sectional profiles", ...
    'subtitle', append(string(size(profiles, 2)), " profiles, looking upstream"), ...
    'min_width', o.min_width, 'max_width', o.max_width)

% (c) elevation of trough and edges along the channel
d = T.dist_along_m/1000;        % [km]
x_lim = [0, max([d; 1e-3])];
ax = nexttile(tl, 3);
hold(ax, 'on')
h_c = plot(ax, d, T.elev_trough_m, 'Color', S.cent, 'LineWidth', S.lw_data);
h_l = plot(ax, d, T.elev_ledge_m, 'Color', S.edge_left, 'LineWidth', S.lw_data);
h_r = plot(ax, d, T.elev_redge_m, 'Color', S.edge_right, 'LineWidth', S.lw_data);
style_axis(ax)
xlim(ax, x_lim)
xlabel(ax, 'Distance along channel [km]', 'Color', S.ink)
ylabel(ax, 'Elevation [m]', 'Color', S.ink)
legend(ax, [h_c, h_l, h_r], {'Centerline (trough)', 'Left edge', 'Right edge'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', S.font_small, 'TextColor', S.ink)
panel_title(ax, "(c)  Elevation along channel", ...
    append("mean depth below left edge: ", sprintf('%.1f', mean(T.depth_m, 'omitnan')), " m"))

% (d) channel width (top) and depth below left edge (bottom, y-axis pointing
% down: deeper channel = lower line), stacked rather than on a right-hand
% y-axis, whose tick labels MATLAB R2026a garbles when printed at high dpi
tl_d = tiledlayout(tl, 2, 1, 'TileSpacing', 'tight', 'Padding', 'none');
tl_d.Layout.Tile = 4;
ax_w = nexttile(tl_d);
plot(ax_w, d, T.width_m/1000, 'Color', S.ink, 'LineWidth', S.lw_data)
style_axis(ax_w)
xlim(ax_w, x_lim)
xticklabels(ax_w, {})
ylabel(ax_w, 'Width [km]', 'Color', S.ink)
panel_title(ax_w, "(d)  Channel width and depth along channel", ...
    append("mean width: ", sprintf('%.2f', mean(T.width_m, 'omitnan')/1000), " km, ", ...
           "mean depth below left edge: ", sprintf('%.1f', mean(T.depth_m, 'omitnan')), " m"))
ax_d = nexttile(tl_d);
plot(ax_d, d, T.depth_m, 'Color', S.ink, 'LineWidth', S.lw_data)
style_axis(ax_d)
ax_d.YDir = 'reverse';
xlim(ax_d, x_lim)
xlabel(ax_d, 'Distance along channel [km]', 'Color', S.ink)
ylabel(ax_d, {'Depth below', 'left edge [m]'}, 'Color', S.ink)
linkaxes([ax_w, ax_d], 'x')

end
