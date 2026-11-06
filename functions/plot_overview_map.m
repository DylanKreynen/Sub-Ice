function plot_overview_map(ax, DEM, R, x_cent, y_cent, edge_coord, x_prof, y_prof, labels, varargin)
%plot_overview_map(ax, DEM, R, x_cent, y_cent, edge_coord, x_prof, y_prof, labels, ...)
%Plots mapped channels on the DEM (see plot_map_base()), in the shared
%Sub-Ice figure style (see fig_style()): centerlines, the channel area
%between the edges (shaded, with outlines along the edges) and, optionally,
%the cross sectional profile transects. The map extent covers all profiles,
%plus a margin.
%
% required input:
% ax = axes handle
% DEM = elevation data array [m] (no data values set to NaN)
% R = spatial referencing information for the DEM [-]
% x_cent, y_cent = cells with centerline coordinates [pix], one per channel
% edge_coord = cell with channel edge coordinates [pix], one per channel (see find_edges())
% x_prof, y_prof = cells with profile coordinates [pix], one per channel (see find_profiles())
% labels = channel labels (string array), one per channel
%
% optional input:
% clims = [min max] elevation for the DEM colors [m] (default: [], automatic)
% hillshade = add a hillshade to the DEM (default: true)
% transects = plot the profile transects (default: false)
% bridge_gaps = connect the edges across gaps left by edge outlier filtering (default: true)
% colors = one color per channel (RGB rows): centerlines in these colors,
%          labelled in the legend instead of on the map, and channel areas
%          not shaded (e.g. for one channel on different dates) (default:
%          [], i.e. black centerlines, shaded channel areas, labels on the map)
% title, subtitle = panel title and subtitle (default: none)
%
% (c) Dylan Kreynen
% University of Oslo
% 2026


%% inputParser

p = inputParser;
addParameter(p, 'clims', [])
addParameter(p, 'hillshade', true)
addParameter(p, 'transects', false)
addParameter(p, 'bridge_gaps', true)
addParameter(p, 'colors', [])
addParameter(p, 'title', "")
addParameter(p, 'subtitle', "")
parse(p, varargin{:});
o = p.Results;


%% actual function

S = fig_style();
no_channels = numel(x_cent);
labels = string(labels);

% map extent: all profiles (first and last point), plus a margin
all_x = cell2mat(cellfun(@(x) reshape(x([1 end], :), [], 1), x_prof(:), 'UniformOutput', false));
all_y = cell2mat(cellfun(@(y) reshape(y([1 end], :), [], 1), y_prof(:), 'UniformOutput', false));
all_x = [all_x; vertcat(x_cent{:})];
all_y = [all_y; vertcat(y_cent{:})];
margin = 0.05*max(max(all_x) - min(all_x), max(all_y) - min(all_y));
plot_map_base(ax, DEM, R, 'clims', o.clims, 'hillshade', o.hillshade, ...
    'extent', [min(all_x)-margin, max(all_x)+margin, min(all_y)-margin, max(all_y)+margin]);

per_channel_colors = ~isempty(o.colors);
h_leg = gobjects(0);
leg_txt = strings(0);

% profile transects (thin line from first to last profile point)
if o.transects
    for c = 1:no_channels
        plot(ax, x_prof{c}([1 end], :), y_prof{c}([1 end], :), 'Color', [S.ink 0.25], 'LineWidth', 0.5, 'HandleVisibility', 'off')
    end
end

if ~per_channel_colors
    % channel areas between the edges, black centerlines, labels at the upstream end
    for c = 1:no_channels
        shade_channel(ax, edge_coord{c}, S.shade, S.shade_alpha, S.shade_outline*S.shade, o.bridge_gaps)
    end
    for c = 1:no_channels
        plot(ax, x_cent{c}, y_cent{c}, 'Color', S.cent, 'LineWidth', S.lw_data, 'HandleVisibility', 'off')
    end
    for c = 1:no_channels
        text(ax, x_cent{c}(1), y_cent{c}(1), append(" ", labels(c), " "), 'Color', S.label, ...
            'FontSize', S.font_small, 'BackgroundColor', 'w', 'Margin', 0.5, 'Interpreter', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
    h_leg = [plot(ax, NaN, NaN, 'Color', S.cent, 'LineWidth', S.lw_data), ...
             patch(ax, NaN, NaN, S.shade, 'FaceAlpha', S.shade_alpha, 'EdgeColor', S.shade_outline*S.shade)];
    leg_txt = ["Centerline", "Channel (between edges)"];

else
    % centerlines and edges (thin) in one color per channel, labelled in the legend
    for c = 1:no_channels
        col = o.colors(c,:);
        for s = [1 3]
            valid = ~isnan(edge_coord{c}(:,s));
            if o.bridge_gaps
                plot(ax, edge_coord{c}(valid,s), edge_coord{c}(valid,s+1), 'Color', [col 0.6], 'LineWidth', 0.6, 'HandleVisibility', 'off')
            else
                plot(ax, edge_coord{c}(:,s), edge_coord{c}(:,s+1), 'Color', [col 0.6], 'LineWidth', 0.6, 'HandleVisibility', 'off')
            end
        end
        h_leg(end+1) = plot(ax, x_cent{c}, y_cent{c}, 'Color', col, 'LineWidth', S.lw_data); %#ok<AGROW>
        leg_txt(end+1) = labels(c); %#ok<AGROW>
    end
    h_leg(end+1) = plot(ax, NaN, NaN, 'Color', [S.muted 0.6], 'LineWidth', 0.6);
    leg_txt(end+1) = "Channel edges";
end

if o.transects
    h_leg(end+1) = plot(ax, NaN, NaN, 'Color', [S.ink 0.25], 'LineWidth', 0.5);
    leg_txt(end+1) = "Cross-sectional profiles";
end
legend(ax, h_leg, leg_txt, 'Location', 'northeast', 'FontSize', S.font_small, 'TextColor', S.ink, 'Box', 'on', 'EdgeColor', S.grid, 'Interpreter', 'none')

if strlength(string(o.title)) > 0
    panel_title(ax, o.title, o.subtitle)
end

end
