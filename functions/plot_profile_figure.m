function plot_profile_figure(ax, prof_dist, profiles, edge_idx, edge_elev, varargin)
%plot_profile_figure(ax, prof_dist, profiles, edge_idx, edge_elev, 'relative', relative, 'title', ttl, ...)
%Plots the cross sectional profiles of one channel, in the shared Sub-Ice
%figure style (see fig_style()): ten profiles along the channel are
%highlighted in color (dark purple upstream to yellow downstream), with
%markers at the detected left and right channel edges; the other profiles
%are shown in light gray in the background. Profiles are shown looking
%upstream, i.e. with the left channel edge (looking downstream) on the
%right hand side.
%
% required input:
% ax = axes handle
% prof_dist = each profile sample's distance from the profile center [m]
%             (positive = left side, looking downstream)
% profiles = matrix containing profiles' sampled elevation [m] (see find_profiles())
% edge_idx = matrix containing profile indices of the channel edges [-] (see find_edges())
% edge_elev = matrix containing channel edge elevations [m] (see find_edges())
%
% optional input:
% relative = plot the profiles between the channel edges only, relative to
%            the elevation of the left channel edge (default: false, i.e.
%            full profiles, absolute elevation)
% title = panel title (default: none)
% subtitle = panel subtitle (default: none)
% min_width, max_width = edge search limits as passed to find_edges() [m],
%            shown as dashed lines (default: 0, i.e. not shown)
% highlight = positions of the highlighted profiles, as normalized distance
%             along the channel [-] (default: 0.05:0.1:0.95, i.e. 10 profiles)
% max_profiles = max. no. of other profiles shown in the background (default: 25)
%
% (c) Dylan Kreynen
% University of Oslo
% 2026


%% inputParser

p = inputParser;
addParameter(p, 'relative', false)
addParameter(p, 'title', "")
addParameter(p, 'subtitle', "")
addParameter(p, 'min_width', 0)
addParameter(p, 'max_width', 0)
addParameter(p, 'highlight', 0.05:0.1:0.95)
addParameter(p, 'max_profiles', 25)
parse(p, varargin{:});

relative = p.Results.relative;
ttl = p.Results.title;
sub = p.Results.subtitle;
min_width = p.Results.min_width;
max_width = p.Results.max_width;
highlight = p.Results.highlight;
max_profiles = p.Results.max_profiles;


%% actual function

S = fig_style();
prof_dist = prof_dist(:);
no_profiles = size(profiles, 2);

% highlighted profiles, and (at most max_profiles) other profiles
hl_idx = unique(min(no_profiles, max(1, round(highlight*no_profiles))));
stride = max(1, ceil(no_profiles/max_profiles));
ctx_idx = setdiff(1:stride:no_profiles, hl_idx);

% edge positions [m] and elevations
edge_pos = NaN(no_profiles, 2);
for s = 1:2
    valid = ~isnan(edge_idx(:,s));
    edge_pos(valid, s) = prof_dist(edge_idx(valid, s));
end
edge_z = edge_elev;

% relative mode: profiles between the edges only, relative to the left edge
if relative
    for i = 1:no_profiles
        if isnan(edge_idx(i,1)) || isnan(edge_idx(i,2))
            profiles(:,i) = NaN;
            continue
        end
        profiles(1:edge_idx(i,2)-1, i) = NaN;
        profiles(edge_idx(i,1)+1:end, i) = NaN;
        profiles(:,i) = profiles(:,i) - edge_elev(i,1);
    end
    edge_z = edge_elev - edge_elev(:,1);
end

hold(ax, 'on')

% centerline and edge search limits (half widths, per side; left = positive)
xline(ax, 0, ':', 'Color', S.muted, 'LineWidth', 0.6, 'HandleVisibility', 'off')
h_lim = [];
for w = {min_width, max_width}
    lim = w{1};
    if isscalar(lim)
        lim = lim/2*[1 1];
    end
    if any(lim > 0)
        h = xline(ax, [lim(1), -lim(2)], '--', 'Color', S.muted, 'LineWidth', 0.6);
        h_lim = h(1);
        set(h(2:end), 'HandleVisibility', 'off')
    end
end

% other profiles (background; none if all profiles are highlighted)
if ~isempty(ctx_idx)
    plot(ax, prof_dist, profiles(:, ctx_idx), 'Color', S.prof_context, 'LineWidth', S.lw_thin, 'HandleVisibility', 'off')
    plot(ax, edge_pos(ctx_idx, :), edge_z(ctx_idx, :), '.', 'Color', S.prof_context*0.85, 'MarkerSize', 6, 'HandleVisibility', 'off')
end

% highlighted profiles
cols = S.cmap_prof_inferno(numel(hl_idx));
h_hl = gobjects(numel(hl_idx), 1);
for k = 1:numel(hl_idx)
    h_hl(k) = plot(ax, prof_dist, profiles(:, hl_idx(k)), 'Color', cols(k,:), 'LineWidth', S.lw_data);
end
h_l = plot(ax, edge_pos(hl_idx, 1), edge_z(hl_idx, 1), 'o', 'MarkerFaceColor', S.edge_left, ...
    'MarkerEdgeColor', S.marker_edge, 'LineWidth', 0.5, 'MarkerSize', 4.5);
h_r = plot(ax, edge_pos(hl_idx, 2), edge_z(hl_idx, 2), 'o', 'MarkerFaceColor', S.edge_right, ...
    'MarkerEdgeColor', S.marker_edge, 'LineWidth', 0.5, 'MarkerSize', 4.5);

% axes
style_axis(ax)
xlim(ax, [min(prof_dist), max(prof_dist)])
xlabel(ax, 'Distance from centerline [m]', 'Color', S.ink)
if relative
    ylabel(ax, 'Elevation relative to left edge [m]', 'Color', S.ink)
else
    ylabel(ax, 'Elevation [m]', 'Color', S.ink)
end
if strlength(string(ttl)) > 0
    panel_title(ax, ttl, sub)
end

% orientation (looking upstream: left edge on the right hand side)
yl = ylim(ax);
y_txt = yl(1) + 0.04*diff(yl);
dx = 0.02*diff(xlim(ax));
text(ax, -dx, y_txt, ['\leftarrow', ' right'], 'Color', S.muted, 'FontSize', S.font_small, 'HorizontalAlignment', 'right')
text(ax, dx, y_txt, ['left ', '\rightarrow'], 'Color', S.muted, 'FontSize', S.font_small, 'HorizontalAlignment', 'left')
text(ax, max(prof_dist) - dx, y_txt, 'view: upstream', 'Color', S.muted, 'FontSize', S.font_small, 'HorizontalAlignment', 'right')

% legend: highlighted profiles (normalized distance along channel), edges, limits
hl_txt = compose("%.2f", hl_idx(:)/no_profiles);
h_leg = [h_hl; h_l; h_r];
leg_txt = [hl_txt; "Left edge"; "Right edge"];
if ~isempty(h_lim)
    h_leg(end+1) = h_lim;
    leg_txt(end+1) = "Min./max. width";
end
lgd = legend(ax, h_leg, leg_txt, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', S.font_small, 'TextColor', S.ink);
title(lgd, {'Norm. distance', 'along channel [-]'}, 'FontSize', S.font_small, 'FontWeight', 'normal', 'Color', S.ink)

end
