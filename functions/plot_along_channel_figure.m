function plot_along_channel_figure(f, T, ttl)
%plot_along_channel_figure(f, T, ttl)
%Plots the channel geometry along one channel, in the shared Sub-Ice
%figure style (see fig_style()): elevation of the trough (at the
%centerline) and of the left and right channel edges (top), channel depth
%below the left edge (middle, y-axis pointing down: deeper channel = lower
%line) and channel width (bottom), against distance along the centerline
%from the channel start point.
%
% required input:
% f = figure handle (the figure is cleared)
% T = table with channel metrics per profile (see channel_metrics())
% ttl = title, e.g. the channel label
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

S = fig_style();
clf(f)
tl = tiledlayout(f, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
d = T.dist_along_m/1000;        % [km]
x_lim = [0, max([d; 1e-3])];

% elevation of trough and edges (two tiles high)
ax1 = nexttile(tl, [2 1]);
hold(ax1, 'on')
h_c = plot(ax1, d, T.elev_trough_m, 'Color', S.cent, 'LineWidth', S.lw_data);
h_l = plot(ax1, d, T.elev_ledge_m, 'Color', S.edge_left, 'LineWidth', S.lw_data);
h_r = plot(ax1, d, T.elev_redge_m, 'Color', S.edge_right, 'LineWidth', S.lw_data);
style_axis(ax1)
xlim(ax1, x_lim)
ylabel(ax1, 'Elevation [m]', 'Color', S.ink)
panel_title(ax1, ttl, append("Mean depth below left edge: ", string(round(mean(T.depth_m, 'omitnan'), 1)), ...
    " m, mean width: ", compose("%.2f", mean(T.width_m, 'omitnan')/1000), " km"))
legend(ax1, [h_c, h_l, h_r], {'Centerline (trough)', 'Left edge', 'Right edge'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', S.font_small, 'TextColor', S.ink)
yl = ylim(ax1);
text(ax1, x_lim(1) + 0.01*diff(x_lim), yl(1) + 0.06*diff(yl), ['downstream ', '\rightarrow'], ...
    'Color', S.muted, 'FontSize', S.font_small)
xticklabels(ax1, {})

% depth below left edge
ax2 = nexttile(tl);
plot(ax2, d, T.depth_m, 'Color', S.ink, 'LineWidth', S.lw_data)
style_axis(ax2)
xlim(ax2, x_lim)
ylabel(ax2, {'Depth below', 'left edge [m]'}, 'Color', S.ink)
ax2.YDir = 'reverse';           % (positive down: deeper channel = lower line)
xticklabels(ax2, {})

% channel width
ax3 = nexttile(tl);
plot(ax3, d, T.width_m/1000, 'Color', S.ink, 'LineWidth', S.lw_data)
style_axis(ax3)
xlim(ax3, x_lim)
ylabel(ax3, 'Width [km]', 'Color', S.ink)
xlabel(ax3, 'Distance along channel [km]', 'Color', S.ink)

linkaxes([ax1, ax2, ax3], 'x')

end
