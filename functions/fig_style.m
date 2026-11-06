function S = fig_style()
%S = fig_style()
%Returns the shared figure style (colors, colormaps, font sizes, line
%widths) used by all Sub-Ice figures, so that all figures look the same.
%Use together with style_axis(), style_colorbar() and panel_title().
%
% Colors have a fixed meaning across all figures:
%   cent                    channel centerline (maps, along-channel plots)
%   shade                   channel area between the edges (maps)
%   edge_left, edge_right   left/right channel edge (looking downstream)
%   prof_highlight          highlighted cross-sectional profiles, upstream to downstream
%   prof_context            other cross-sectional profiles, in the background
%
% output:
% S = struct with colors (RGB, 0-1), colormap functions, font sizes and
%     line widths
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

% neutrals
S.ink    = [ 31  31  30]/255;     % primary ink: lines, text
S.muted  = [107 106  99]/255;     % secondary ink: axes, tick labels, subtitles
S.grid   = [228 226 218]/255;     % grid lines
S.excl   = [236 235 230]/255;     % shaded (excluded) areas

% channels
S.cent       = [0 0 0];           % centerline (black)
S.shade      = [204 121 167]/255; % channel area between the edges, on maps (light purple, transparent)
S.shade_alpha = 0.30;             % opacity of the channel area shading on maps (the DEM shows through)
S.shade_outline = 0.65;           % outline of the channel area shading: shade color times this (darker shade)
S.edge_left  = [106  44 138]/255; % left channel edge (dark purple)
S.edge_right = [190 140 215]/255; % right channel edge (light purple)
S.marker_edge = [0 0 0];          % marker outline (black)

% maps
S.label      = [0 0 0];           % channel labels (black, on a white background)

% colormaps (n = no. of colors)
% DEM on maps: light grays only (no black), so the channel shading stays clear over low areas
S.cmap_dem_light = @(n) repmat(linspace(0.45, 0.97, n)', 1, 3);
% profiles along channel (start to end), e.g. when all profiles are colored
S.cmap_prof = @(n) parula(n);
% highlighted profiles along channel (inferno 0.2 ... 0.9, dark purple upstream to yellow downstream)
inferno_ramp = [66 10 104; 106 23 110; 147 38 103; 188 55 84; 221 81 58; 243 120 25; 252 165 10; 246 215 70]/255;
S.cmap_prof_inferno = @(n) interp1(linspace(0, 1, size(inferno_ramp, 1)), inferno_ramp, linspace(0, 1, n));
S.prof_context = [0.87 0.87 0.87];  % other profiles (and their markers) behind the highlighted ones (light gray)
% channels, e.g. when several channels are plotted together (categorical)
S.cmap_channels = @(n) lines(n);

% fonts [pt]
S.font = 8;                       % default text
S.font_small = 7;                 % legends, annotations
S.font_title = 9;                 % panel titles

% line widths [pt]
S.lw_data = 1.2;                  % data lines
S.lw_thin = 0.8;                  % secondary lines

end
