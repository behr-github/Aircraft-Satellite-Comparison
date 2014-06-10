function [no2,alt,utc,interp_height] = boundary_layer_verification( Merge, Data, tz, varargin )
%[ pix_lat, pix_lon, BEHR_NO2, aircraft_NO2, OMI_NO2 ]

%boundary_layer_verification: Compare satellite and aircraft measurements
%using the boundary layer method. Returns pixel latitude, longitude, BEHR
%tropospheric NO2, aircraft trop. NO2, and OMI trop. NO2.
%   This function compares satellite and aircraft NO2 columns using the
%   boundary layer method described in Russell et. al. Atmos. Chem. Phys.
%   11, 8543-8554, 2011. As inputs, it requires a Merge data structure (the
%   result of reading an ICART file using read_icart_file.m), a BEHR Data
%   data structure (one of the outputs of BEHR_main.m), and a time zone
%   abbreviation:
%
%       EST = Eastern Std.      
%       CST = Central Std.      
%       MST = Mountain Std.     
%       PST = Pacific Std.      
%
%   Since Aura satellite overpass is ~1:45 p in standard time, do not use
%   daylight savings time values.
%
% Parameters
%   timerange: By default, this method will restrict flight data to between
%   12:00p and 3:00p local, as per Russell et. al.  This can be overridden
%   using this parameter, which requires a cell array with the start and
%   end times as strings in military time format.
%
%   NO2_sep: This is the minimum difference between adjacent NO2
%   measurements to qualify as a boundary layer entrance exit.  It defaults
%   to 1000 pptv.
%
%   time_sep: The maximum time (in seconds) between "clusters" of boundary
%   layer entrances and exits - see "utc_sep" in
%   findall_no2_bdy_layer_heights.m for more information. Defaults to 1000 s.
%
%   cloud_product: Which cloud product (omi or modis) to use in rejecting
%   pixels.  Defaults to omi.
%
%   cloud_frac_max: The maximum allowed geometric cloud fraction.  Defaults
%   to 0.2; recommended value for use with MODIS cloud product is 0.
%
%   rowanomaly: The method of rejecting pixels based on row anomaly.
%   Defaults to 'AlwaysByRow'.  See omi_pixel_reject or omi_rowanomaly for
%   more information on the possible choices ('AlwaysByRow', 'RowsByTime',
%   'XTrackFlags', and 'XTrackFlagsLight').

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst',})));
p.addParamValue('timerange',{'12:00','15:00'},@iscell)
p.addParamValue('NO2_sep',1e3,@isscalar);
p.addParamValue('time_sep',1e3,@isscalar);
p.addParamValue('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis'})));
p.addParamValue('cloud_frac_max',0.2, @isscalar);
p.addParamValue('rowanomaly','AlwaysByRow',@(x) strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'}));

p.parse(Merge,Data,tz,varargin{:});
pout = p.Results;

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
timerange = pout.timerange;
NO2_sep = pout.NO2_sep;
time_sep = pout.time_sep;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;

% Load the aircraft data
no2 = Merge.Data.NO2_UCB.Values;
alt = Merge.Data.ALTP.Values;
utc = Merge.Data.UTC.Values;
lat = Merge.Data.LATITUDE.Values;
lon = Merge.Data.LONGITUDE.Values;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  RESTRICT TO TIME RANGE  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove any values outside the time range (12:00-3:00p local standard by
% default)
utcstart = local2utc(timerange{1}, tz); utcend = local2utc(timerange{2}, tz);
time_logical = utc >= utcstart & utc <= utcend;
no2 = no2(time_logical);
alt = alt(time_logical);
utc = utc(time_logical);
lat = lat(time_logical);
lon = lon(time_logical);

% Find any fill, ULOD, or LLOD values in any of the imported data sets and
% remove them.
no2_fill = Merge.Data.NO2_UCB.Fill;
alt_fill = Merge.Data.ALTP.Fill;
lat_fill = Merge.Data.LATITUDE.Fill;
lon_fill = Merge.Data.LONGITUDE.Fill;

fill_logical = no2 ~= no2_fill & alt ~= alt_fill & lat ~= lat_fill & lon ~= lon_fill;

ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;

lod_logical = no2 ~= ulod & no2 ~= llod & alt ~= ulod & alt ~= llod; % Lat, lon, and utc should not be subject to limits of detection

no2 = no2(lod_logical & fill_logical);
alt = alt(lod_logical & fill_logical);
utc = utc(lod_logical & fill_logical);
lat = lat(lod_logical & fill_logical);
lon = lon(lod_logical & fill_logical);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  FIND BDY LAYER HEIGHTS  %%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[heights, times] = findall_no2_bdy_layer_heights(no2, NO2_sep, alt, utc, time_sep, 'NO2');
if isempty(heights)
    error('bdy_layer_verify:findBoundaryLayer','Could not find boundary layer heights.')
end

interp_height = interp1(times, heights, utc); % Linearly interpolate the boundary layer height to every value of UTC.

% For values outside of the range of "times," assume that the boundary
% layer height equals the closest value.
first_height = find(~isnan(interp_height),1,'first');
last_height = find(~isnan(interp_height),1,'last');

if first_height > 1;
    interp_height(1:first_height-1) = interp_height(first_height);
end
if last_height < numel(interp_height)
    interp_height(last_height+1:end) = interp_height(last_height);
end





end

