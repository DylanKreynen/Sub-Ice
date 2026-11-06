function DEM_smooth = smooth_dem(DEM, res, window_DEM, path_to_DEM)
%DEM_smooth = smooth_dem(DEM, res, window_DEM, path_to_DEM)
%Returns a smoothed version of the DEM (2D moving mean), to be used for
%centerline mapping. The smoothed DEM is cached to disk (next to the input
%DEM) and loaded from there in subsequent runs with the same window size.
%
% required input:
% DEM = elevation data array [m] (no data values set to NaN)
% res = resolution of DEM [m/pix]
% window_DEM = window size for DEM smoothing [m] (will be rounded up to
%              [pix], set to 0 for no smoothing)
% path_to_DEM = path to the DEM GeoTIFF, used to construct the filename of
%               the cached smoothed DEM ("<DEM name>_smooth_<window_DEM>.mat")
%
% output:
% DEM_smooth = smoothed elevation data array [m] (same size as DEM; equal
%              to DEM if window_DEM = 0)
%
% (c) Dylan Kreynen
% University of Oslo
% 2024-2026


if window_DEM == 0
    % no smoothing requested
    DEM_smooth = DEM;
    return
end

window = ceil(window_DEM/res);  % from [m] to [pix]
[DEM_dir, DEM_name] = fileparts(path_to_DEM);
fn_smooth = fullfile(DEM_dir, append(DEM_name, '_smooth_', string(window_DEM), '.mat'));

if isfile(fn_smooth)
    % load cached smoothed DEM
    DEM_smooth = load(fn_smooth).DEM_smooth;
    disp(append("Loaded smoothed version of DEM from disk: ", fn_smooth, ". Window size: ", string(window_DEM), " [m] (ca. ", string(window), " [pix]). "))

    % check if DEM dimensions agree (e.g. cache from a different DEM)
    if ~isequal(size(DEM), size(DEM_smooth))
        error(append("Dimensions of DEM and cached smoothed version (", fn_smooth, ") do not agree. Delete the cached file to smooth the DEM again. "))
    end

else
    % smooth DEM and cache to disk
    disp(append("Could not find smoothed version of DEM. Smoothing now. Window size: ", string(window_DEM), " [m] (ca. ", string(window), " [pix]). "))
    DEM_smooth = smoothdata2(DEM, 'movmean', window, 'omitnan');
    save(fn_smooth, 'DEM_smooth', '-v7.3');     % -v7.3 required for variables > 2 GB
    disp(append("Saved smoothed version of DEM to disk: ", fn_smooth, ". "))
end

end
