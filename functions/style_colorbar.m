function style_colorbar(hcb, label)
%style_colorbar(hcb, label)
%Applies the shared Sub-Ice figure style (see fig_style()) to a colorbar.
%
% required input:
% hcb = colorbar handle
%
% optional input:
% label = colorbar label, shown along the bar
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

S = fig_style();
hcb.Color = S.muted;
hcb.TickDirection = 'out';
hcb.LineWidth = 0.6;
hcb.FontSize = S.font;
hcb.Box = 'off';
if nargin > 1
    ylabel(hcb, label, 'Color', S.muted)
end

end
