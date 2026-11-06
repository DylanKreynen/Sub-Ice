function [DEM, R] = read_dem(path_to_DEM, DEM_nodata)
%[DEM, R] = read_dem(path_to_DEM, DEM_nodata)
%Reads a DEM from a GeoTIFF and sets no data values to NaN. The no data
%value is read from the GeoTIFF (if given there); DEM_nodata (e.g. from the
%configuration file) is used as well, as fallback.
%
% required input:
% path_to_DEM = path to the DEM (GeoTIFF, projected coordinates in meters)
%
% optional input:
% DEM_nodata = no data value of the DEM (default: [], i.e. only the value
%              given in the GeoTIFF)
%
% output:
% DEM = elevation data array [m] (single or double precision, no data = NaN)
% R = spatial referencing information for the DEM [-]
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

if nargin < 2
    DEM_nodata = [];
end

[DEM, R] = readgeoraster(path_to_DEM);

% integer DEMs: to floating point, so no data can be set to NaN
if isinteger(DEM)
    DEM = single(DEM);
end

% no data values: as given in the GeoTIFF, and DEM_nodata
nodata = georasterinfo(path_to_DEM).MissingDataIndicator;
nodata = unique([nodata(:); DEM_nodata(:)]);
nodata = nodata(~isnan(nodata));
for k = 1:numel(nodata)
    DEM(DEM == nodata(k)) = NaN;
end

end
