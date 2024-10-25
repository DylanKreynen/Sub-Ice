function [x_samp, y_samp] = centerline_query_pts(P_start, P_end, res, search_step, search_angle, no_samp_pts, first_after_start)
%[x_samp, y_samp] = centerline_query_pts(P_start, P_end, res, search_step, search_angle, no_samp_pts, first_after_start)
%Alternative function for "centerline_query_pts_line()" that works in a 
%similar way, but returns coordinates of sampling points on a circular arc 
%transect rather than a straight line transect. 
%
% input: 
% P_start = vector containing x and y coordinates of start point [pix coord]
% P_end = vector containing x and y coordinates of end point [pix coord]
% res = image resolution [m/pixel]
% search_step = distance to step away from P_start to construct search profile [m]
% search_angle = angle of view within to look for centerline [deg]
% no_samp_pts = number of sampling points on search profile [-] (will be made odd number)
% first_after_start = 1 or 0, whether we're looking for the first point 
% after channel start point
%
% output: 
% x_samp = vector with x coordinates of sampling locations
% y_samp = vector with y coordinates of sampling locations
% 
% Dylan Kreynen (University of Oslo)
% June - Oct 2024

x_start = P_start(1); 
y_start = P_start(2); 
x_end = P_end(1); 
y_end = P_end(2); 

search_step = search_step/res;          % from [m] to [pix]
search_angle = deg2rad(search_angle);   % [rad]

search_dir = [x_end-x_start, y_end-y_start];    % perpendicular direction to start-end vector (search direction)
search_dir = search_dir./norm(search_dir);      % now a unit vector

search_start = atan(search_dir(2)/search_dir(1)) + search_angle/2;   % in [rad]
search_stop = atan(search_dir(2)/search_dir(1)) - search_angle/2;    % in [rad]

no_samp_pts = ceil(no_samp_pts/2);                              % to end up with uneven number of samp pts (see next line)
th = linspace(search_start, search_stop, 2*no_samp_pts+1);      % always one point directly in search direction

if first_after_start == 1
    if search_dir(1)>=0 && search_dir(2)>0 || search_dir(1)>=0 && search_dir(2)<0   % in quadrants I or IV
        x_samp = search_step*cos(th) + x_start;
        y_samp = search_step*sin(th) + y_start;
    else                                                                            % in quadrants II or III
        x_samp = -search_step*cos(th) + x_start;
        y_samp = -search_step*sin(th) + y_start;
    end
else
    if search_dir(1)>=0 && search_dir(2)>0 || search_dir(1)>=0 && search_dir(2)<0   % in quadrants I or IV
        x_samp = search_step*cos(th) + x_end;
        y_samp = search_step*sin(th) + y_end;
    else                                                                            % in quadrants II or III
        x_samp = -search_step*cos(th) + x_end;
        y_samp = -search_step*sin(th) + y_end;
    end
end

end