function [res_x, idx_of_result] = knee_pt(y, varargin)
% [res_x, idx_of_result] = knee_pt(y, x, just_return, knee_method)
% Returns the x-location of a (single) knee of curve y=f(x) (this is 
% useful for e.g. figuring out where the eigenvalues peter out), and 
% returns the index of the x-coordinate at the knee. 
%
% required inputs:
% y = vector of the y values, at least three elements long
%
% optional inputs:
% x = vector of x values the same size as y (default: incremental array of
% same size as y)
% knee_method = method used to detect knee point,  "LinearRegression" or "Kneedle" (default: "LinearRegression")
%
%
% Important:  
%  The x and y  don't need to be sorted, they just have to
%  correspond: knee_pt([1,2,3],[3,4,5]) = knee_pt([3,1,2],[5,3,4])
%
%  Because of the way the function operates y must be at least 3
%  elements long and the function will never return either the first or the
%  last point as the answer.
%
%
% The "Kneedle" method operates using a Kneedle algorithm, which draws a line
% between the start and end points, then draws perpendicular lines from
% that line to fit to each point. it then selects the knee point as the
% maximum length line. 
%
% The "LinearRegression" method operates by walking along the curve one bisection point at a time and
% fitting two lines, one to all the points to left of the bisection point and one
% to all the points to the right of of the bisection point.
% The knee is judged to be at a bisection point which minimizes the
% sum of errors for the two fits.
%
% The errors being used are sum(abs(del_y)) or RMS depending on the
% (very obvious) internal switch.  Experiment with it if the point returned
% is not to your liking -- it gets pretty subjective...
%
%
% 
% Original function by Dmitry Kaplan
% retrieved from MATLAB Central File Exchange
% 2012
% 
% Kneedle method by Adam Isaak 
% University of Manitoba
% 2026

%% input parser
% default parameter values
default_x = NaN;
default_knee_method = "LinearRegression";

% parse input arguments
p = inputParser;
validXY = @(x) isvector(x) && length(x) >= 3 && isnumeric(x);
validKneeMethod = @(x) convertCharsToStrings(x)=="LinearRegression" | convertCharsToStrings(x)=="Kneedle";
addRequired(p, 'y', validXY);
addOptional(p, 'x', default_x, validXY);
addOptional(p, 'knee_method', default_knee_method, validKneeMethod);
parse(p, y, varargin{:});

y = p.Results.y;
x = p.Results.x;
knee_method = p.Results.knee_method;

% Verify that both x and y are the same dimension
if isnan(x)
    x = 1:length(y);
end

x = x(:);
y = y(:);

if ~all(size(x) == size(y))
    error('knee_pt: y and x must have the same dimensions');
end

%% actual function

%set internal operation flags
use_absolute_dev_p = true;  %ow quadratic


% default return values
res_x = nan;
idx_of_result = nan;


% make sure the x and y are sorted in increasing X-order
if any(diff(x) < 0)
    [~,idx]=sort(x);
    y = y(idx);
    x = x(idx);
else
    idx = 1:length(x);
end


if knee_method == "LinearRegression"
    %the code below "unwraps" the repeated regress(y,x) calls.  It's
    %significantly faster than the former for longer y's
    %
    %figure out the m and b (in the y=mx+b sense) for the "left-of-knee"
    sigma_xy = cumsum(x.*y);
    sigma_x  = cumsum(x);
    sigma_y  = cumsum(y);
    sigma_xx = cumsum(x.*x);
    n        = (1:length(y))';
    det = n.*sigma_xx-sigma_x.*sigma_x;
    mfwd = (n.*sigma_xy-sigma_x.*sigma_y)./det;
    bfwd = -(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det;
    %figure out the m and b (in the y=mx+b sense) for the "right-of-knee"
    sigma_xy = cumsum(x(end:-1:1).*y(end:-1:1));
    sigma_x  = cumsum(x(end:-1:1));
    sigma_y  = cumsum(y(end:-1:1));
    sigma_xx = cumsum(x(end:-1:1).*x(end:-1:1));
    n        = (1:length(y))';
    det = n.*sigma_xx-sigma_x.*sigma_x;
    mbck = flipud((n.*sigma_xy-sigma_x.*sigma_y)./det);
    bbck = flipud(-(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det);
    %figure out the sum of per-point errors for left- and right- of-knee fits
    error_curve = nan(size(y));
    for breakpt = 2:length(y-1)
        delsfwd = (mfwd(breakpt).*x(1:breakpt)+bfwd(breakpt))-y(1:breakpt);
        delsbck = (mbck(breakpt).*x(breakpt:end)+bbck(breakpt))-y(breakpt:end);
        %disp([sum(abs(delsfwd))/length(delsfwd), sum(abs(delsbck))/length(delsbck)])
        if (use_absolute_dev_p)
            % error_curve(breakpt) = sum(abs(delsfwd))/sqrt(length(delsfwd)) + sum(abs(delsbck))/sqrt(length(delsbck));
            error_curve(breakpt) = sum(abs(delsfwd))+ sum(abs(delsbck));
        else
            error_curve(breakpt) = sqrt(sum(delsfwd.*delsfwd)) + sqrt(sum(delsbck.*delsbck));
        end
    end
    %find location of the min of the error curve
    [~,loc] = min(error_curve);
    res_x = x(loc);
    idx_of_result = idx(loc);
    
elseif knee_method == "Kneedle"
    line_dist = nan(size(y));

    % find the line in the y=mx+b
    start_pos_x = x(1);
    start_pos_y = y(1);
    end_pos_x = x(end);
    end_pos_y = y(end);
    m = (end_pos_y - start_pos_y)/(end_pos_x - start_pos_x);
    b = start_pos_y - m*(start_pos_x);

    % calculate the angle between the line and each cross section point
    theta = (pi()/2) - atan(m);
    
    % for each calculate the length of each perpendicular line from
    % start-end line back to the current cross section point
    for i = 2:(length(y)-1)
        tri_pos = (m*x(i))+b;
        h = y(i) - tri_pos;
        line_dist(i) = sin(theta)*h;
    end

    % return index of the point with the longest perpendicular line 
    [~, loc] = max(line_dist);
    res_x = x(loc);
    idx_of_result = idx(loc);
end

end
