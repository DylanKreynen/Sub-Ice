function f = new_figure(size_cm, hidden)
%f = new_figure(size_cm, hidden)
%Creates a figure for the extended Sub-Ice figures (one or more per
%channel), so that a run does not open a large number of figure windows:
%figures that are only saved to file are not shown (close them after
%saving), figures that are not saved are shown docked, as tabs in a single
%figure window.
%
% required input:
% size_cm = figure size [width height] [cm] (used when saved; [] = MATLAB default)
% hidden = true: do not show the figure (e.g. when it is saved to file);
%          false: show the figure, docked
%
% output:
% f = figure handle
%
% (c) Dylan Kreynen
% University of Oslo
% 2026

if hidden
    f = figure('Visible', 'off');
else
    f = figure('WindowStyle', 'docked');
end
if ~isempty(size_cm)
    f.Units = 'centimeters';
    f.Position(3:4) = size_cm;
end

end
