function shade_channel(ax, edge_coord, col, alpha, outline_col, bridge_gaps)
%shade_channel(ax, edge_coord, col, alpha, outline_col, bridge_gaps)
%Shades the channel area between the left and right channel edges on a map,
%with thin outlines along the left and right edges (not across the channel
%ends).
%
% required input:
% ax = axes handle (map, in image coordinates [pix])
% edge_coord = matrix containing channel edge coordinates [pix] (see find_edges()),
%              one row per profile: [x_left y_left x_right y_right]
% col = fill color (RGB)
% alpha = fill opacity [0-1]
% outline_col = color of the edge outlines (RGB)
%
% optional input:
% bridge_gaps = connect the edges across profiles where an edge was not
%               detected or left out as outlier (default: false, i.e. gaps in
%               outlines and shading)
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

if nargin < 6
    bridge_gaps = false;
end

valid_l = ~isnan(edge_coord(:,1));
valid_r = ~isnan(edge_coord(:,3));
valid = valid_l & valid_r;

if bridge_gaps
    % one polygon over all profiles where both edges are available
    i = find(valid);
    if numel(i) >= 2
        patch(ax, [edge_coord(i,1); flipud(edge_coord(i,3))], [edge_coord(i,2); flipud(edge_coord(i,4))], col, ...
            'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    end
    plot(ax, edge_coord(valid_l,1), edge_coord(valid_l,2), 'Color', outline_col, 'LineWidth', 0.6, 'HandleVisibility', 'off')
    plot(ax, edge_coord(valid_r,3), edge_coord(valid_r,4), 'Color', outline_col, 'LineWidth', 0.6, 'HandleVisibility', 'off')

else
    % one polygon per run of profiles where both edges are available
    d = diff([0; valid; 0]);
    run_start = find(d == 1);
    run_end = find(d == -1) - 1;
    for r = 1:numel(run_start)
        i = run_start(r):run_end(r);
        if numel(i) < 2
            continue
        end
        patch(ax, [edge_coord(i,1); flipud(edge_coord(i,3))], [edge_coord(i,2); flipud(edge_coord(i,4))], col, ...
            'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    end
    % (NaN = not available: gap in the line)
    plot(ax, edge_coord(:,1), edge_coord(:,2), 'Color', outline_col, 'LineWidth', 0.6, 'HandleVisibility', 'off')
    plot(ax, edge_coord(:,3), edge_coord(:,4), 'Color', outline_col, 'LineWidth', 0.6, 'HandleVisibility', 'off')
end

end
