function [] = lines_to_shp(x, y, R, path, varargin)
%[] = lines_to_shp(x, y, R, path, field, value)
%Takes a lines (x and y in image coordinates [pix]) and writes it to a
%georeferenced shapefile for use in e.g. QGIS or ArcGIS. 
%
% required input: 
% x = matrix containing lines' x coordinates as columns [pix]
% y = matrix containing lines' y coordinates as columns [pix]
% R = spatial referencing information for the image array [-] (from e.g. readgeoraster)
% path = path where to save file (incl. filename but without extension) [-]
% 
% optional input: 
% field = string, header for column in shapefile's attribute table (dynamic property, one per shapefile)
% value = vector with values to save in said attribute table column (one per line)
% attributes = table with additional attribute table columns, one row per line
%              (in the order the lines are written; numeric, string, char or
%              cellstr columns; names should be max. 10 characters for .dbf)
%
% note: 'x', 'y', and 'value' are also allowed to be of type 'cell' of equal length, 
% in which case they may contain multiple matrics (but must all be of type 'cell')
%
% output: 
% shapefile (written to disk)
% 
% (c) Dylan Kreynen
% University of Oslo
% 2024-2026

%% inputParser

% default parameter values
default_field = 'line_no'; 
default_value =  1:size(x, 1);  

% parse input arguments
p = inputParser; 
validMapCellsRef = @(x) class(x) == "map.rasterref.MapCellsReference";  
validStringChar = @(x) isstring(x) | ischar(x); 
addRequired(p, 'x')
addRequired(p, 'y')
addRequired(p, 'R', validMapCellsRef)
addRequired(p, 'path', validStringChar)
addOptional(p, 'field', default_field, validStringChar)
addOptional(p, 'value', default_value)
addParameter(p, 'attributes', table(), @istable)
parse(p, x, y, R, path, varargin{:});

field = p.Results.field;
value = p.Results.value;
attributes = p.Results.attributes;


%% actual function

if ~iscell(x)
    % if not a cell, make it a cell
    x = {x}; 
    y = {y}; 
    value = {value}; 
end

[x_world, y_world] = intrinsicToWorld(R, x{1}(:,1), y{1}(:,1)); % image to map coordinates
gshape = geoshape(y_world, x_world, field, value{1}(1)); 

if size(x{1}, 2) > 1
    % add columns (lines) as additional lines to geoshape
    for i = 2:size(x{1}, 2)
        [x_world, y_world] = intrinsicToWorld(R, x{1}(:,i), y{1}(:,i)); 
        gshape = append(gshape, y_world, x_world, field, value{1}(i)); 
    end
end

if length(x) > 1
    % add cells (channels) as additional lines to shapefile
    for c = 2:length(x)
        [x_world, y_world] = intrinsicToWorld(R, x{c}(:,1), y{c}(:,1)); % image to map coordinates
        gshape = append(gshape, y_world, x_world, field, value{c}(1)); 

        if size(x{c}, 2) > 1
            % add columns (lines) as additional lines to geoshape
            for i = 2:size(x{c}, 2)
                [x_world, y_world] = intrinsicToWorld(R, x{c}(:,i), y{c}(:,i));
                gshape = append(gshape, y_world, x_world, field, value{c}(i)); 
            end  
        end
    end
end

% add additional attributes (one row per line)
if ~isempty(attributes)
    if height(attributes) ~= length(gshape)
        error(append("Number of attribute rows (", string(height(attributes)), ") does not match number of lines (", string(length(gshape)), "). "))
    end
    names = attributes.Properties.VariableNames;
    for a = 1:length(names)
        val = attributes.(names{a});
        if isstring(val) || ischar(val)
            val = cellstr(val);
        end
        gshape.(names{a}) = val';
    end
end

% .shp, .shx, .dbf files:
dbfspec = makedbfspec(gshape);
shapewrite(gshape, path, 'DbfSpec', dbfspec);

% .prj file: 
wkt = wktstring(R.ProjectedCRS); 
writematrix(wkt, append(path, ".prj"), 'FileType', 'text', 'QuoteStrings', false); 