function [ line_x, line_y, legend_str, P,R ] = calc_fit_line( x, y, varargin )
%plot_fit_line Draws a fit line and 1:1 line on the current figure.
%   Calculates a line of best fit for the given data after removing NaNs.
%   Returns vectors of x & y points to plot the line, a string formatted
%   with information on the fit line for a legend, the polynomial
%   coefficients ( [slope, intercept] ) and the R^2 value.
%
%   By default this uses the normal polyfit function to determine the fit
%   line, to use a reduced major axis regression (function by Edward
%   Peitzer, http://www.mbari.org/staff/etp3/regress/index.htm), set the
%   parameter 'regression' to 'majoraxis'.  This is good when both x & y
%   have errors, rather than assuming the x values are 100% correct. Set
%   the parameter 'regression' to 'rma' to use his geometric mean fit
%   (lsqfitgma).

p = inputParser;
p.addRequired('x',@isnumeric);
p.addRequired('y',@isnumeric);
p.addParameter('one2one',true);
p.addParameter('regression','y-resid',@(x) any(strcmpi(x,{'y-resid','majoraxis','RMA'})));

p.parse(x,y,varargin{:});
pout = p.Results;
x = pout.x;
y = pout.y;
one2one = pout.one2one;
regression = pout.regression;

if any(isnan(x)) || any(isnan(y))
    warning('NaNs detected, removing any points with a value of NaN for either coordinate');
    nans = isnan(x) | isnan(y);
    x = x(~nans); y = y(~nans);
end

if strcmpi(regression,'y-resid');
    [P, R] = polyfit_R2(x,y,1);
elseif strcmpi(regression,'majoraxis');
    [P(1), P(2), R] = lsqfitma(x,y);
    R = R^2; % lsqfitma does not square R before returning.
elseif strcmpi(regression,'RMA');
    [P(1), P(2), R] = lsqfitgm(x,y);
    R = R^2;
end

line_x = 0:1e16:1e17;
line_y = polyval(P,line_x);
legend_str = sprintf('Fit: %.4fx + %.2g \nR^2 = %.4f',P(1),P(2),R);

end

