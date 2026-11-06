function panel_title(ax, ttl, sub)
%panel_title(ax, ttl, sub)
%Left-aligned bold panel title, e.g. '(a)  M2', with an optional muted
%subtitle, in the shared Sub-Ice figure style (see fig_style()).
%
% required input:
% ax = axes handle
% ttl = title
%
% optional input:
% sub = subtitle
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

S = fig_style();
title(ax, ttl, 'Color', S.ink, 'FontSize', S.font_title, 'FontWeight', 'bold')
if nargin > 2 && any(strlength(string(sub)) > 0)     % (sub can be multi-line: cell array)
    subtitle(ax, sub, 'Color', S.muted, 'FontSize', S.font_small + 0.5)
end
ax.TitleHorizontalAlignment = 'left';

end
