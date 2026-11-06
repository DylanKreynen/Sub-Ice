function style_axis(ax, grid_dir)
%style_axis(ax, grid_dir)
%Applies the shared Sub-Ice figure style (see fig_style()) to an axes: no
%box, outward ticks, muted axis color, light grid lines behind the data,
%left-aligned title.
%
% required input:
% ax = axes handle
%
% optional input:
% grid_dir = grid lines to draw: 'y' (default), 'x', 'xy' or 'none'
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

if nargin < 2
    grid_dir = 'y';
end
S = fig_style();

box(ax, 'off')
if strcmp(grid_dir, 'none')
    ax.Layer = 'top';           % (no grid: axis lines on top, e.g. not covered by an image's edge)
else
    ax.Layer = 'bottom';        % (grid lines behind the data)
end
ax.TickDir = 'out';
ax.LineWidth = 0.6;
ax.FontSize = S.font;
ax.XColor = S.muted;
ax.YColor = S.muted;
ax.XGrid = contains(grid_dir, 'x');
ax.YGrid = contains(grid_dir, 'y');
ax.GridColor = S.grid;
ax.GridAlpha = 1;
ax.GridLineStyle = '-';
ax.Title.Color = S.ink;
ax.Title.FontSize = S.font_title;
ax.TitleHorizontalAlignment = 'left';

end
