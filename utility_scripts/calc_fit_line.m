function [ line_x, line_y, legend_str, LineData ] = calc_fit_line( x, y, varargin )
%plot_fit_line Draws a fit line and 1:1 line on the current figure.
%   Calculates a line of best fit for the given data after removing NaNs.
%   Returns vectors of x & y points to plot the line, a string formatted
%   with information on the fit line for a legend, and a structure
%   containing the slope/intercept vector P (compatible with polyval), the
%   R^2 values, and std. deviations of slope and intercept.
%
%   This makes use of the various linear regression functions by Edward
%   Peitzer (http://www.mbari.org/staff/etp3/regress/index.htm). Which
%   regression is used can be set by the 'regression' parameter.
%       
%       'y-resid' - (default) minimizes differences between the points and
%       the fit line in the y-direction only.  Best when x is considered to
%       have no error.
%
%       'x-resid' - minimized differences between the points and the fit
%       line in the x-direction only.
%
%       'majoraxis' - minimized the normal distance between the fit line
%       and the points. Said to fit a line along the major axis of an
%       ellipse that represents the variance.  Should only be used if x and
%       y have the same units and scale, and if both x & y have
%       non-negligible error.
%
%       'RMA' - takes a geometric mean of the y-on-x and x-on-y regression.
%       Effectively minimizes the area of a triangle between each point and
%       the fit line. Most computationally expensive, but good when x & y
%       may not have the same scale or units and both have error terms.
%
%   http://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov5.html
%   has a good explanation of the various sorts of regressions.

p = inputParser;
p.addRequired('x',@isnumeric);
p.addRequired('y',@isnumeric);
p.addParameter('regression','y-resid',@(x) any(strcmpi(x,{'y-resid','x-resid','majoraxis','RMA'})));

p.parse(x,y,varargin{:});
pout = p.Results;
x = pout.x;
y = pout.y;
regression = lower(pout.regression); % explicitly make the regression string lower case to ease comparison in the switch-case statement

if any(isnan(x)) || any(isnan(y))
    warning('NaNs detected, removing any points with a value of NaN for either coordinate');
    nans = isnan(x) | isnan(y);
    x = x(~nans); y = y(~nans);
end

sem = -1e9;
switch regression
    case 'y-resid'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfity(x,y);
        R = R^2; % All of the lsqfit** functions do not square R before returning
        
    case 'x-resid'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfitx(x,y);
        
    case 'majoraxis'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfitma(x,y);
        R = R^2; 
        
    case 'rma'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfitgm(x,y);
        R = R^2;
end

LineData.P = P;
LineData.R2 = R;
LineData.StdDevM = sigma_m;
LineData.StdDevB = sigma_b;

line_x = 0:1e16:1e17;
line_y = polyval(P,line_x);
legend_str = sprintf('Fit: %.4fx + %.2g \nR^2 = %.4f',P(1),P(2),R);

end

