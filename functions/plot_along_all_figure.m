function plot_along_all_figure(f, T_along, mean_elev, labels, varargin)
%plot_along_all_figure(f, T_along, mean_elev, labels, 'colors', colors, 'title', ttl, 'subtitle', sub)
%Plots mean profile elevation, channel depth below the left edge (y-axis
%pointing down: deeper channel = lower line) and channel width against
%distance along the centerline, for several channels (or one channel on
%several dates) together, in the shared Sub-Ice figure style (see
%fig_style()).
%
% required input:
% f = figure handle (the figure is cleared)
% T_along = cell with tables with channel metrics per profile, one per
%           channel (see channel_metrics())
% mean_elev = cell with mean elevation of each profile [m], one per channel
% labels = legend labels (string array), one per channel
%
% optional input:
% colors = one color per channel (RGB rows) (default: categorical colors)
% title, subtitle = title and subtitle of the top panel (default: none)
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

S = fig_style();
no_channels = numel(T_along);

p = inputParser;
addParameter(p, 'colors', S.cmap_channels(no_channels))
addParameter(p, 'title', "")
addParameter(p, 'subtitle', "")
parse(p, varargin{:});
cols = p.Results.colors;

clf(f)
tl = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(3, 1);
y_lab = {'Mean profile elevation [m]', 'Depth below left edge [m]', 'Width [km]'};
for k = 1:3
    ax(k) = nexttile(tl);
    hold(ax(k), 'on')
    for c = 1:no_channels
        d = T_along{c}.dist_along_m/1000;
        switch k
            case 1, y = mean_elev{c}(:);
            case 2, y = T_along{c}.depth_m;
            case 3, y = T_along{c}.width_m/1000;
        end
        plot(ax(k), d, y, 'Color', cols(c,:), 'LineWidth', S.lw_data)
    end
    style_axis(ax(k))
    ylabel(ax(k), y_lab{k}, 'Color', S.ink)
    if k == 2
        ax(k).YDir = 'reverse';     % (depth positive down: deeper channel = lower line)
    end
    if k < 3
        xticklabels(ax(k), {})
    end
end
xlabel(ax(3), 'Distance along channel [km]', 'Color', S.ink)
linkaxes(ax, 'x')
xlim(ax(1), [0, max(cellfun(@(T) max([T.dist_along_m; 1]), T_along))/1000])
if strlength(string(p.Results.title)) > 0
    panel_title(ax(1), p.Results.title, p.Results.subtitle)
end
xl = xlim(ax(1));
yl = ylim(ax(1));
text(ax(1), xl(1) + 0.01*diff(xl), yl(1) + 0.08*diff(yl), ['downstream ', '\rightarrow'], ...
    'Color', S.muted, 'FontSize', S.font_small, 'VerticalAlignment', 'middle')
legend(ax(1), string(labels), 'Location', 'northoutside', 'Orientation', 'horizontal', ...
    'Box', 'off', 'FontSize', S.font_small, 'TextColor', S.ink, 'NumColumns', min(no_channels, 6), 'Interpreter', 'none')

end
