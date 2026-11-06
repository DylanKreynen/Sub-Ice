function T_sum = channel_summary(T, label, end_found, channel_length)
%T_sum = channel_summary(T, label, end_found, channel_length)
%Returns a one-row table with summary statistics of a channel, from its
%channel metrics per profile (see channel_metrics()).
%
% required input:
% T = table with channel metrics per profile (see channel_metrics())
% label = channel label (e.g. channel name, or date for a time series)
% end_found = 1 if the centerline connects start and end point [-]
% channel_length = centerline length [m]
%
% output:
% T_sum = table (one row) with columns: channel, end_found, length_km,
%         no_profiles, mean/median/min/max/std width [m], mean left and
%         right half width [m], mean/median/max depth below left edge [m],
%         mean left edge, right edge and trough elevation [m], and the
%         percentage of left edges, right edges and troughs flagged as
%         outlier [%]
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

T_sum = table(string(label), end_found, channel_length/1000, height(T), ...
    mean(T.width_m, 'omitnan'), median(T.width_m, 'omitnan'), min(T.width_m), max(T.width_m), std(T.width_m, 'omitnan'), ...
    mean(T.hw_left_m, 'omitnan'), mean(T.hw_right_m, 'omitnan'), ...
    mean(T.depth_m, 'omitnan'), median(T.depth_m, 'omitnan'), max(T.depth_m), ...
    mean(T.elev_ledge_m, 'omitnan'), mean(T.elev_redge_m, 'omitnan'), mean(T.elev_trough_m, 'omitnan'), ...
    100*mean(~T.valid_ledge), 100*mean(~T.valid_redge), 100*mean(~T.valid_trough), ...
    'VariableNames', {'channel', 'end_found', 'length_km', 'no_profiles', ...
    'mean_width_m', 'median_width_m', 'min_width_m', 'max_width_m', 'std_width_m', ...
    'mean_hw_left_m', 'mean_hw_right_m', ...
    'mean_depth_m', 'median_depth_m', 'max_depth_m', ...
    'mean_elev_ledge_m', 'mean_elev_redge_m', 'mean_elev_trough_m', ...
    'outliers_ledge_pct', 'outliers_redge_pct', 'outliers_trough_pct'});

end
