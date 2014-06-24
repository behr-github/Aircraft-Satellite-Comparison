function [ bin_values, bin_midpoints, bin_error ] = bin_rolling_vertical_profile( altitude, data_values, binwidth, binspacing, varargin )
%bin_rolling_vertical_profile(altitude, data, binwidth, binspacing, [options]) Bin vertical profile data with rolling median or mean
%   This function groups data into altitude bins of width "binwidth" (in
%   km) and with centers "binspacing" km apart, essentially a rolling
%   median or rolling mean. By default, this function finds the bin value as the
%   median, with error presented as 25th and 75th quantiles.  To use mean
%   and std. deviation instead, enter 'mean' as the optional fourth
%   argument.
%
%   This function will return 3 variables: bin_values, the data values
%   median or mean for the bin; bin_midpoints, the midpoint altitude of the
%   bin; and bin_error, a 2 x n matrix of quartiles if running in median
%   mode, or a 1 x n row vector of std. dev. in mean mode.
%
%   Josh Laughner <joshlaugh5@gmail.com> 24 June 2014

p = inputParser;
p.addRequired('altitude',@isnumeric);
p.addRequired('data_values',@isnumeric);
p.addRequired('binwidth',@isscalar);
p.addRequired('binspacing',@isscalar);
p.addOptional('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));

p.parse(altitude, data_values, binwidth, binspacing, varargin{:});
pout = p.Results;
altitude = pout.altitude;
data_vals = pout.data_values;
binwidth = pout.binwidth;
binspacing = pout.binspacing;
binmode = pout.binmode;

% Define the edges and centers of the bins.  The first bin center will be
% 1/2 binwidth above 0.

top_alt = max(altitude(:));
bin_tops = binwidth:binspacing:100; % Generate more bins than we could possibly need...
xx = bin_tops <= (top_alt + binwidth); %...then restrict then to those that are needed to hold our data
bins = [bin_tops(xx) - binwidth; bin_tops(xx) - 0.5*binwidth; bin_tops(xx)]'; % Create a matrix where the columns are [bin_bottom, bin_center, bin_top]

% Do the binning. 
bin_values = zeros(size(bins,1),1);
if strcmpi(binmode,'mean')
    bin_error = zeros(size(bins,1),1);
else
    bin_error = zeros(size(bins,1),2);
end

for a=1:numel(bin_values)
    bin_data_vals = data_vals(altitude > bins(a,1) & altitude <= bins(a,3));
    if strcmpi(binmode,'mean')
        bin_values(a) = nanmean(bin_data_vals(:));
        bin_error(a) = nanstd(bin_data_vals(:)) / sqrt(numel(bin_data_vals));
    else
        bin_values(a) = nanmedian(bin_data_vals(:));
        bin_error(a,:) = quantile(bin_data_vals(:),[0.25,0.75]);
    end
end

bin_midpoints = bins(:,2);

end

