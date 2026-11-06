function keep = check_points(R, xs, ys, xe, ye, labels, varargin)
%keep = check_points(R, xs, ys, xe, ye, labels, 'DEM', DEM, 'shp', path_to_shp)
%Checks channel start and end points against the DEM, before mapping:
% - points outside the DEM: the channel is skipped (with a warning); if no
%   channel is left, an error is given (most likely, the points are in a
%   different projection than the DEM)
% - points on DEM no data: warning (the centerline search will likely fail)
% - shapefile projection (.prj) different from the DEM's: warning
%
% required input:
% R = spatial referencing information for the DEM [-]
% xs, ys, xe, ye = start and end points, image coordinates of the DEM [pix]
%                  (one per channel)
% labels = channel labels (or dates), one per channel
%
% optional input:
% DEM = elevation data array [m] (for the no data check; default: [], no check)
% shp = path to the shapefile the points (or centerlines) were read from
%       (for the projection check; default: '', no check)
% what = what is checked, for the messages (default: "channel")
%
% output:
% keep = logical, one per channel: true if start and end point lie within the DEM
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

p = inputParser;
addParameter(p, 'DEM', [])
addParameter(p, 'shp', '')
addParameter(p, 'what', "channel")
parse(p, varargin{:});
DEM = p.Results.DEM;
shp = p.Results.shp;
what = string(p.Results.what);
labels = string(labels(:));

% projection of the shapefile vs. DEM
if strlength(string(shp)) > 0 && isfile(shp) && isprop(R, 'ProjectedCRS') && ~isempty(R.ProjectedCRS)
    crs_shp = shapeinfo(shp).CoordinateReferenceSystem;
    if ~isempty(crs_shp) && ~isequal(crs_shp, R.ProjectedCRS)
        disp(append("Warning: the projection of ", shp, " (", crs_shp.Name, ") seems to differ from the DEM's (", ...
            R.ProjectedCRS.Name, "). Points and lines should be in the same projection as the DEM. "))
    end
end

% points within the DEM extent
inside = @(x, y) x >= 0.5 & x <= R.RasterSize(2) + 0.5 & y >= 0.5 & y <= R.RasterSize(1) + 0.5;
keep = inside(xs(:), ys(:)) & inside(xe(:), ye(:));
if ~any(keep)
    error(append("All ", what, " start and end points lie outside the DEM. Check that they are in the same projection as the DEM (map coordinates in meters). "))
elseif any(~keep)
    disp(append("Warning: start or end point outside the DEM for ", what, " ", strjoin(labels(~keep), ", "), ...
        ": skipped. Check its start and end point (same projection as the DEM?). "))
end

% points on DEM no data
if ~isempty(DEM)
    z = @(x, y) DEM(sub2ind(size(DEM), min(max(round(y), 1), size(DEM, 1)), min(max(round(x), 1), size(DEM, 2))));
    on_nodata = false(size(keep));
    on_nodata(keep) = isnan(z(xs(keep), ys(keep))) | isnan(z(xe(keep), ye(keep)));
    if any(on_nodata)
        disp(append("Warning: start or end point on DEM no data for ", what, " ", strjoin(labels(on_nodata), ", "), ...
            ". The centerline search will likely fail; consider moving the point(s) onto the DEM. "))
    end
end

end
