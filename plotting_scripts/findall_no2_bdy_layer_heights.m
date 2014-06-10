function [ heights, median_times ] = findall_no2_bdy_layer_heights( data, data_sep, altitude, utc, utc_sep, field_name)
%findall_bdy_layer_heights: Output boundary layer heights for a full day's
%worth of data.
%   This function will analyze a full day of aircraft data and find the
%   median boundary height for clusters of boundary layer crossings.  This
%   is based on Ashley Russell's method for comparing satellite data
%   against boundary layer flights, described in her 2011 paper.
%
%   This function operates by identifying clusters of boundary layer
%   crossings that occur in a relatively short period of time. (The utc_sep
%   input determines the maximum allowed time between crossings that belong
%   to the same cluster). The subset of data defined by this cluster is
%   then passed to the "find_bdy_layer_heights" function, which calulates
%   the boundary layer height for that cluster based on its criteria.
%
%   This function requires 6 inputs:
%       data: the measurements of the parameter for which you are finding
%       the boundary layer, e.g. NO2 or potential temperature.
%
%       data_sep: The difference in value between two successive
%       measurements that corresponds to a possible boundary layer entrance
%       or exit. For NO2, ~1000 is a good number to start with.
%
%       altitude: the altitude of the plane in kilometers.  This must have
%       a 1-to-1 correspondence with the data input.
%
%       utc: The times that correspond to both the altitude and the data
%       input.  This can be any time data, although UTC is perhaps the most
%       consistent.
%
%       utc_sep: The time gap that the function uses to separate clusters
%       of boundary layer crossings.  1000 s is usually a good starting
%       guess.
%
%       field_name: The name of the Merge structure field that "data" was
%       read from.  This is necessary as the function
%       "find_bdy_layer_height" (called below) uses difference criteria for
%       different fields.
%
%   The values returned are the calculated median boundary layer heights
%   and the times corresponding to those heights.

yy = find(abs(diff(data)) > data_sep);

% If we are looking at an NO2 field, reject the differences if
%   a) The concurrent NO2 measurement is greater than 10,000 (this usually
%   means the aircraft is passing through a plume from a city, forest fire,
%   etc. that might confuse the boundary layer location)
%   b) There is no consistent change in altitude as well
if any(strcmp(field_name,{'NO2_LIF','NO2_UCB','NO2','NO2_MixingRatio_LIF'}));
%     % Check if either value that contributed to the difference is greater
%     % than 10,000 pptv
%     tmp_no2 = [data(yy); data(yy+1)];
%     tmp = tmp_no2 < 1e4;
%     yy_high_no2 = false(1,size(tmp,2));
%     for a = 1:size(tmp,2)
%         yy_high_no2(a) = tmp(1,a) & tmp(2,a);
%     end
%     yy = yy(yy_high_no2);
    
    % Check that the difference in altitude that corresponds to the change
    % in NO2 is greater than 20 m/s
    dz_dt = diff(data)./diff(utc);
    yy = yy(dz_dt(yy) > 0.02);
end

zz = find(abs(diff(utc(yy+1))) > utc_sep);
n = numel(zz) - 1;

heights = -9e9*ones(1,n);
median_times = -9e9*ones(1,n);

for a=1:n
    %if any(strcmp(field_name,{'NO2_LIF','NO2_UCB','NO2'})) && any(data(yy(zz(a))+1:yy(zz(a+1)))>=1e4)
    %else
    blh = find_bdy_layer_height(data(yy(zz(a))+1:yy(zz(a+1))),altitude(yy(zz(a))+1:yy(zz(a+1))), field_name);
    t = median(utc(yy(zz(a))+1:yy(zz(a+1))));
    
    heights(a) = blh;
    median_times(a) = t;
    %end
end

heights(heights < -1e9) = []; median_times(median_times < -1e9) = [];

end

