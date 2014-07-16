function [ P,R ] = plot_fit_line( x, y, varargin )
%plot_fit_line Draws a fit line and 1:1 line on the current figure.
%   Automatically draws a fit line on the current plot (given the data on
%   the plot as input variables - i.e. after plot(x,y) call
%   plot_fit_line(x,y) to make the lines). Also creates a legend with the
%   linear fit info.  By default also draws a 1:1 line.  This can be
%   disabled by setting the parameter 'one2one' to 0.
%
%   By default this uses the normal polyfit function to determine the fit
%   line, to use a reduced major axis regression (function by Edward
%   Peitzer, http://www.mbari.org/staff/etp3/regress/index.htm), set the
%   parameter 'regression' to 'majoraxis'.  This is good when both x & y
%   have errors, rather than assuming the x values are 100% correct.

p = inputParser;
p.addRequired('x',@isnumeric);
p.addRequired('y',@isnumeric);
p.addParamValue('one2one',true);
p.addParamValue('regression','lineary',@(x) any(strcmpi(x,{'lineary','majoraxis'})));

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

if strcmpi(regression,'lineary');
    [P, R] = polyfit_R2(x,y,1);
elseif strcmpi(regression,'majoraxis');
    [P(1), P(2), R] = lsqfitma(x,y);
    R = R^2; % lsqfitma does not square R before returning.
end
line(0:1e16:1e17,polyval(P,0:1e16:1e17),'color','k','linestyle','--','linewidth',2);
legendcell = {'All points',sprintf('Fit: %.4fx + %.2g \nR^2 = %.4f',P(1),P(2),R)};

if one2one
    line(0:1e16:1e17,1:1e16:1e17,'color','r','linestyle',':','linewidth',2);
    legendcell{end+1} = '1:1';
end
legend(legendcell{:})
end

