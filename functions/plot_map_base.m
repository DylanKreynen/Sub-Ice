function clims = plot_map_base(ax, DEM, R, varargin)
%clims = plot_map_base(ax, DEM, R, 'clims', clims, 'hillshade', hillshade, 'extent', extent)
%Draws the DEM as background of a map (in light grays, optionally with a
%hillshade), with an elevation colorbar, a scalebar and a north arrow, in
%the shared Sub-Ice figure style (see fig_style()). The map is drawn in
%image coordinates of the (full resolution) DEM [pix], so that centerlines,
%edges etc. can be plotted on top directly. Large DEMs are subsampled for
%display.
%
% required input:
% ax = axes handle
% DEM = elevation data array [m] (no data values set to NaN)
% R = spatial referencing information for the DEM [-]
%
% optional input:
% clims = [min max] elevation for the colormap [m] (default: [], i.e. 2nd
%         and 98th percentile of the elevations within the map extent)
% hillshade = add a hillshade to the DEM (default: true)
% hs_strength = strength of the hillshade [0-1] (default: 0.7)
% hs_z = vertical exaggeration of the hillshade [-] (default: 40; relief on
%        ice shelves is only a few m)
% extent = [x_min x_max y_min y_max] map extent in image coordinates [pix]
%          (default: [], i.e. the full DEM)
% max_pix = max. no. of DEM pixels to display along each side; larger DEMs
%           are subsampled for display (default: 1500)
%
% output:
% clims = [min max] elevation of the colormap [m]
%
% (c) Dylan Kreynen
% University of Oslo
% 2026


%% inputParser

p = inputParser;
addParameter(p, 'clims', [])
addParameter(p, 'hillshade', true)
addParameter(p, 'hs_strength', 0.7)
addParameter(p, 'hs_z', 40)
addParameter(p, 'extent', [])
addParameter(p, 'max_pix', 1500)
parse(p, varargin{:});

clims = p.Results.clims;
hillshade_on = p.Results.hillshade;
hs_strength = p.Results.hs_strength;
hs_z = p.Results.hs_z;
extent = p.Results.extent;
max_pix = p.Results.max_pix;


%% actual function

S = fig_style();
res = R.CellExtentInWorldX;         % resolution of DEM [m/pix]

% map extent [pix], within the DEM
if isempty(extent)
    extent = [1 size(DEM, 2) 1 size(DEM, 1)];
end
c0 = max(1, floor(extent(1)));  c1 = min(size(DEM, 2), ceil(extent(2)));
r0 = max(1, floor(extent(3)));  r1 = min(size(DEM, 1), ceil(extent(4)));

% subsample for display
f = max(1, ceil(max(c1-c0+1, r1-r0+1)/max_pix));
D = double(DEM(r0:f:r1, c0:f:c1));
x_ext = [c0, c0 + (size(D, 2)-1)*f];      % pixel centers, full resolution image coordinates
y_ext = [r0, r0 + (size(D, 1)-1)*f];

% color limits
if isempty(clims)
    v = sort(D(~isnan(D)));                 % (2nd and 98th percentile, without toolbox)
    clims = v(max(1, round([0.02 0.98]*numel(v))))';
    if isempty(clims) || clims(1) == clims(2)
        clims = [min(D(:)), max(D(:)) + 1];
    end
end

% DEM in light grays, optionally times a hillshade, as an RGB image (the
% colormap set below still gives the elevation colorbar)
cmap = S.cmap_dem_light(256);
g = interp1(linspace(clims(1), clims(2), 256), cmap(:,1), min(max(D, clims(1)), clims(2)));
if hillshade_on && hs_strength > 0
    % (faded out over 10 m above the colorbar's top, so e.g. rough grounded ice stays plain)
    strength = hs_strength*min(max((clims(2) + 10 - D)/10, 0), 1);
    g = min(1, g.*(1 + strength.*(hillshade(D, hs_z, res*f) - 1)));
end
g(isnan(g)) = 1;
image(ax, x_ext, y_ext, repmat(g, 1, 1, 3), 'AlphaData', ~isnan(D))
ax.Color = 'w';                             % (no data: white)
hold(ax, 'on')
axis(ax, 'image')
xlim(ax, x_ext + [-0.5 0.5]*f)
ylim(ax, y_ext + [-0.5 0.5]*f)
colormap(ax, cmap)
clim(ax, clims)
set(ax, 'XTick', [], 'YTick', [], 'XColor', S.muted, 'YColor', S.muted, 'FontSize', S.font, 'Layer', 'top')
box(ax, 'on')
hcb = colorbar(ax);
style_colorbar(hcb, 'Elevation [m]')

% scalebar (bottom left, ca. 1/5 of the map width) and north arrow above it
xl = xlim(ax);
yl = ylim(ax);
nice = [1 2 5]'*10.^(-1:3);
nice = sort(nice(:));
sb_km = nice(find(nice*1000/res <= diff(xl)/5, 1, 'last'));
if isempty(sb_km)
    sb_km = nice(1);
end
sb_pix = sb_km*1000/res;
x0 = xl(1) + 0.05*diff(xl);
y0 = yl(2) - 0.05*diff(yl);
plot(ax, [x0, x0+sb_pix], [y0, y0], 'w', 'LineWidth', 5, 'HandleVisibility', 'off')
plot(ax, [x0, x0+sb_pix], [y0, y0], 'Color', S.ink, 'LineWidth', 2, 'HandleVisibility', 'off')
text(ax, x0+sb_pix/2, y0 - 0.015*diff(yl), append(string(sb_km), " km"), 'Color', S.ink, ...
    'BackgroundColor', 'w', 'Margin', 1, 'FontSize', S.font_small, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

% north arrow: true north at the center of the map (needs the DEM's
% projected coordinate reference system)
if isprop(R, 'ProjectedCRS') && ~isempty(R.ProjectedCRS)
    [xw, yw] = intrinsicToWorld(R, mean(xl), mean(yl));
    [lat, lon] = projinv(R.ProjectedCRS, xw, yw);
    [xn, yn] = projfwd(R.ProjectedCRS, min(lat + 0.01, 89.99), lon);
    [xi, yi] = worldToIntrinsic(R, xn, yn);
    u = [xi - mean(xl), yi - mean(yl)];
    if norm(u) > 0
        u = u/norm(u);                          % (image coordinates, y down)
        len = 0.06*diff(yl);
        base =[x0 + sb_pix/2, y0 - 0.08*diff(yl) - 0.5*len] - 0.5*len*u;
        tip = base + len*u;
        n = [-u(2), u(1)];
        head = [tip; tip - 0.35*len*u + 0.18*len*n; tip - 0.35*len*u - 0.18*len*n];
        plot(ax, [base(1) tip(1)], [base(2) tip(2)], 'Color', S.ink, 'LineWidth', 1.5, 'HandleVisibility', 'off')
        patch(ax, head(:,1), head(:,2), S.ink, 'EdgeColor', S.ink, 'HandleVisibility', 'off')
        text(ax, tip(1) + 0.3*len*u(1), tip(2) + 0.3*len*u(2), 'N', 'FontWeight', 'bold', 'Color', S.ink, ...
            'FontSize', S.font, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    end
end

end


function shade = hillshade(Z, z_factor, cell_size)
% relative hillshade of DEM Z (y down): 1 on flat ground, < 1 facing away from
% and > 1 facing the light (from the top left of the map, 45 deg above the horizon)
Z = fillmissing(Z, 'nearest');                              % (no gaps along the ice front and no-data areas)
[gx, gy] = gradient(z_factor*Z, cell_size);                 % (gy along the rows, i.e. down)
L = [-sind(45)*cosd(45), cosd(45)*cosd(45), sind(45)];      % light direction (x right, y up, z up)
shade = (-gx*L(1) + gy*L(2) + L(3))./sqrt(gx.^2 + gy.^2 + 1)/L(3);
shade = min(max(shade, 0), 1.5);
end
