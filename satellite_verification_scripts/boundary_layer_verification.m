function [omi_lon, omi_lat, OMI_NO2, BEHR_NO2, DC8_NO2] = boundary_layer_verification( Merge, Data, tz, varargin )
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
%   daylight savings time values.  Further, you must pass only one element
%   of Data; the one corresponding to the swath of interest.
%
% Parameters
%   timerange: By default, this method will restrict flight data to between
%   12:00p and 3:00p local, as per Russell et. al.  This can be overridden
%   using this parameter, which requires a cell array with the start and
%   end times as strings in military time format.
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
%
%   DEBUG_LEVEL: The level of output messages to write; 0 = none, 1 =
%   normal; 2 = verbose

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst',})));
p.addParamValue('timerange',{'12:00','15:00'},@iscell)
p.addParamValue('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis'})));
p.addParamValue('cloud_frac_max',0.2, @isscalar);
p.addParamValue('rowanomaly','AlwaysByRow',@(x) strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'}));
p.addParamValue('DEBUG_LEVEL',1,@isscalar);

p.parse(Merge,Data,tz,varargin{:});
pout = p.Results;

% Check that only one element of Data was passed
if numel(Data)>1; error('bdy_layer_verify:DataInput','Only pass one top-level element of the Data structure'); end

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
timerange = pout.timerange;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

% Load the aircraft data
if DEBUG_LEVEL > 0; fprintf(' Getting data and restricting to 12:00-15:00 %s\n',tz); end
no2_raw = Merge.Data.NO2_UCB.Values;
alt_raw = Merge.Data.ALTP.Values;
radar_alt_raw = Merge.Data.Radar_Altitude.Values;
utc_raw = Merge.Data.UTC.Values;
lat_raw = Merge.Data.LATITUDE.Values;
lon_raw = Merge.Data.LONGITUDE.Values - 360;
pressure_raw = Merge.Data.PRESSURE.Values;
temperature_raw = Merge.Data.TEMPERATURE.Values;

% Get the date
merge_date = Merge.metadata.date;

% Load the file containing the altitude ranges in which to look for
% boundary layer crossings. Loads the variable "Ranges," a structure.

load('/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/ARCTAS-CA Altitude Ranges Exclusive.mat');
% Store the date corresponding to each index for easy searching later
range_dates = {datestr(Ranges(1).Date,29), datestr(Ranges(2).Date,29), datestr(Ranges(3).Date,29), datestr(Ranges(4).Date,29)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  RESTRICT TO TIME RANGE  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove any values outside the time range (12:00-3:00p local standard by
% default)
utcstart = local2utc(timerange{1}, tz); utcend = local2utc(timerange{2}, tz);
time_logical = utc_raw >= utcstart & utc_raw <= utcend;
no2 = no2_raw(time_logical);
alt = alt_raw(time_logical);
radar_alt = radar_alt_raw(time_logical);
utc = utc_raw(time_logical);
lat = lat_raw(time_logical);
lon = lon_raw(time_logical);
pressure = pressure_raw(time_logical);

% We do not want to restrict temperature, as it will be used to get a
% temperature profile for conversion of 

% Find any fill, ULOD, or LLOD values in any of the imported data sets relating to NO2
% and remove them.
no2_fill = Merge.Data.NO2_UCB.Fill;
alt_fill = Merge.Data.ALTP.Fill;
ralt_fill = Merge.Data.Radar_Altitude.Fill;
lat_fill = Merge.Data.LATITUDE.Fill;
lon_fill = Merge.Data.LONGITUDE.Fill;
pres_fill = Merge.Data.PRESSURE.Fill;

fill_logical = no2 ~= no2_fill & alt ~= alt_fill & lat ~= lat_fill & lon ~= lon_fill & pressure ~= pres_fill & radar_alt == ralt_fill;

ulod = Merge.metadata.upper_lod_flag;
llod = Merge.metadata.lower_lod_flag;

lod_logical = no2 ~= ulod & no2 ~= llod & alt ~= ulod & alt ~= llod & pressure ~= ulod & pressure ~= llod; % Lat, lon, radar altitude, and utc should not be subject to limits of detection

no2 = no2(lod_logical & fill_logical);
alt = alt(lod_logical & fill_logical);
radar_alt = radar_alt(lod_logical & fill_logical);
utc = utc(lod_logical & fill_logical);
lat = lat(lod_logical & fill_logical);
lon = lon(lod_logical & fill_logical);
pressure = pressure(lod_logical & fill_logical);

% Now for temperature
temperature_fill = Merge.Data.TEMPERATURE.Fill;
T_fill_logical = temperature_raw ~= temperature_fill & alt_raw ~= alt_fill;
T_lod_logical = temperature_raw ~= ulod & temperature_raw ~= llod & alt_raw ~= ulod & alt_raw ~= llod;
temperature = temperature_raw(T_fill_logical & T_lod_logical);
alt_T = alt_raw(T_fill_logical & T_lod_logical);
utc_T = utc_raw(T_fill_logical & T_lod_logical);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  FIND BDY LAYER HEIGHTS  %%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DEBUG_LEVEL > 0; fprintf(' Calculating & interpolating boundary layer heights.\n'); end

% Get the range for today;
d = find(datenum(range_dates) == datenum(merge_date));

[heights, times] = findall_no2_bdy_layer_heights(utc_raw, no2_raw, alt_raw, Ranges(d).Ranges);
nans = isnan(heights);
heights = heights(~nans); times = times(~nans);
if numel(heights) < 2
    error('bdy_layer_verify:findBoundaryLayer','Could not find 2 boundary layer heights - needed for interpolation')
end

interp_height = interp1(times, heights, utc_raw); % Linearly interpolate the boundary layer height to every value of UTC.

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

% Now clip the interpolated height to the time range and clear values
% corresponding to NO2 fills/lod flags
interp_height = interp_height(time_logical);
interp_height = interp_height(lod_logical & fill_logical);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  CALCULATE NO2 COLUMNS   %%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integrate aircraft NO2 measurements from the boundary layer height to the
% surface (i.e. over the distance defined by radar altitude).  If the plane
% is above the boundary layer, set column to NaN since no data is available
% on BL NO2.
if DEBUG_LEVEL > 0; fprintf(' Integrating columns\n'); end

% First, extract a temperature profile; we will need this to convert mixing
% ratios to number densities
[T_bins, T_bin_alt] = bin_vertical_profile(alt_T, temperature, 0.3);
% The temperature profiles often have an inversion in the BL, so split the
% temperature data into slope <0 (remnant PBL/free troposphere) and >=0
% (mixed layer)
ML_ind = find(diff(T_bins) < 0,1,'first')+1; % The +1 accounts for the fact that length(diff(x)) is length(x) - 1
[T_poly_mixed, T_R2_mixed] = polyfit_R2(T_bin_alt(1:ML_ind),T_bins(1:ML_ind),1);
[T_poly_free, T_R2_free] = polyfit_R2(T_bin_alt(ML_ind:end),T_bins(ML_ind:end),1);
% We'll calculate the mixed layer height as the intersection of the two
% fits
A = [-T_poly_free(1), 1; -T_poly_mixed(1), 1]; b = [T_poly_free(2); T_poly_mixed(2)];
X = A\b; ML_height = X(1)*1000;

% We need the P0 and scale height (H) used in the conversion of pressure to
% altitude.  Fit the data to a line such that ln(P) = (-1/H)z + ln(P0),
% thus H = -1/slope and P0 = e^(intercept)
P = polyfit(alt_raw,log(pressure_raw),1);
P0 = exp(P(2)); H = -1/P(1);

% Also, define Avogadro's number and gas constant;
Av = 6.022e23; % molec. / mol
R = 8.314e4; % (hPa * cm^3) / (mol * K)

% Get the GLOBE terrain altitude for this region
latminmax = [floor(min(lat)), ceil(max(lat))];
lonminmax = [floor(min(lon)), ceil(max(lon))];
globe_dir = '/Volumes/share/GROUP/SAT/BEHR/GLOBE_files';
[terpres, refvec] = globedem(globe_dir,1,latminmax,lonminmax);
cell_count = refvec(1);
globe_latmax = refvec(2); globe_latmin = globe_latmax - size(terpres,1)*(1/cell_count);
globe_lat_matrix = (globe_latmin + 1/(2*cell_count)):(1/cell_count):globe_latmax;
globe_lat_matrix = globe_lat_matrix';
globe_lat_matrix = repmat(globe_lat_matrix,1,size(terpres,2));

globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(terpres,2)*(1/cell_count);
globe_lon_matrix = globe_lonmin + 1/(2*cell_count):(1/cell_count):globe_lonmax;
globe_lon_matrix = repmat(globe_lon_matrix,size(terpres,1),1); 

terpres(isnan(terpres)) = 0; % The Globe database fills oceans with NANs.  Since we care about column height, we want oceans to be at 0 altitude.

% Prepare variables from the BEHR Data Structure
omi_lat = Data.Latitude; omi_lon = Data.Longitude;
TP_pres = Data.TropopausePressure;

% Now numerically integrate the columns.
no2_integrated_col = -127*ones(1,numel(no2));
dz = 1; % width of integration step in meters
for c = 1:numel(no2)
    
    if alt(c) > interp_height(c);
        if DEBUG_LEVEL > 1; fprintf('\t Column %d of %d - above BL\n',c,numel(no2)); end
        no2_column = NaN; %Skip any columns where the plane is above the BL
    elseif no2(c) < 0;
        no2_column = NaN; %Also skip any columns where the NO2 measurement is negative
    else
        if DEBUG_LEVEL > 1; fprintf('\t Column %d of %d\n',c,numel(no2)); end
        % We'll integrate in groups of 1 m (100 cm)
        no2_column = 0;
        BLH = round(interp_height(c)*1000); %Get the boundary layer height in meters
        
        % Find the OMI pixel that corresponds to this measurement, get the
        % tropopause altitude
        dlat = abs(omi_lat - lat(c)); dlon = abs(omi_lon - lon(c));
        dtot = dlat + dlon; xx = (dtot == min(dtot(:)));
        TP_alt = round(log(TP_pres(xx)/P0)*(-H)*1000); % Converted to meters
        
        %Radar altitude is the distance from the plane to surface, pressure
        %altitude is the planes "absolute" altitude, therefore the surface
        %altitude should be the difference.  Convert to meters.
%        surface_alt = round((alt(c) - radar_alt(c))*1000); 
        
        % Since the radar was apparently broken on all four ARCTAS flights,
        % we'll instead find the closest GLOBE measurement and use that for
        % surface altitude
        
        dlat = abs(globe_lat_matrix - lat(c)); dlon = abs(globe_lon_matrix - lon(c));
        dtot = dlat + dlon; xx = (dtot == min(dtot(:)));
        surface_alt = round(terpres(xx)); % Globe measurements should be in meters above sea level
        
        if surface_alt < BLH
            % Surface to boundary layer
            for h=surface_alt:dz:BLH
                P = P0 * exp(-h/(H*1000));
                if h <= ML_height; % Account for the temperature inversion in the mixed layer
                    %T = polyval(T_poly_mixed,h/1000);
                    T = T_poly_mixed(1)*(h/1000) + T_poly_mixed(2);
                elseif h > ML_height;
                    %T = polyval(T_poly_free,h/1000);
                    T = T_poly_free(1)*(h/1000) + T_poly_free(2);
                end
                conc_NO2 = (Av * P * no2(c)*1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
            % Boundary layer to tropopause
            for h=(BLH+1):dz:TP_alt
                P = P0 * exp(-h/(H*1000));
                %T = polyval(T_poly_free,h/1000);
                T = T_poly_free(1)*(h/1000) + T_poly_free(2);
                conc_NO2 = (Av * P * 40*1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv.  Assume 40 pptv in the free troposphere
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
        else
            % Surface to tropopause
            for h=surface_alt:dz:TP_alt
                P = P0 * exp(-h/(H*1000));
                %T = polyval(T_poly_free,h/1000);
                T = T_poly_free(1)*(h/1000) + T_poly_free(2);
                conc_NO2 = (Av * P * 40*1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv.  Assume 40 pptv in the free troposphere
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
        end
    end
    no2_integrated_col(c) = no2_column;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% AVERAGE COLUMNS TO PIXELS  %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DEBUG_LEVEL > 0; fprintf(' Averaging columns to pixels\n'); end

% Find pixels that fail the criteria used for BEHR mapping.
% omi_pixel_reject.m works by setting Areaweight to 0 for any pixels that
% fail the criteria. Hence, we need to create an areaweight field that
% assumes all pixels are valid.

Data.Areaweight = ones(size(Data.Longitude));
if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
Data2 = omi_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
xx = Data2.Areaweight > 0;

% Read in corner lat/lon and BEHR NO2 columns. Latcorn/loncorn are 4 x n
% fields, where the first dimension is each corner for a pixel.
omi_lat = Data2.Latitude(xx); omi_lon = Data2.Longitude(xx);
corner_lat = Data2.Latcorn(:,xx); corner_lon = Data2.Loncorn(:,xx);
BEHR_NO2 = Data2.BEHRColumnAmountNO2Trop(xx);
OMI_NO2 = Data2.ColumnAmountNO2Trop(xx);

% Remove from consideration all pixels whose corners fall entirely outside
% the range of lat/lons that the aircraft flew
if DEBUG_LEVEL > 0; fprintf('   Removing pixels outside aircraft track\n'); end
lon_logical = corner_lon > min(lon) & corner_lon < max(lon);
lat_logical = corner_lat > min(lat) & corner_lat < max(lat);
latlon_logical = any(lon_logical) & any(lat_logical); % any() operates on matrices along the first dimension; here that is each corner for pixel n

omi_lat = omi_lat(latlon_logical);
omi_lon = omi_lon(latlon_logical);
corner_lat = corner_lat(:,latlon_logical);
corner_lon = corner_lon(:,latlon_logical);
BEHR_NO2 = BEHR_NO2(latlon_logical);
OMI_NO2 = OMI_NO2(latlon_logical);

% For each pixel that remains, average the aircraft columns to each pixel
if DEBUG_LEVEL > 0; fprintf('   Averaging columns\n'); end
DC8_NO2 = NaN(size(BEHR_NO2));
for a=1:numel(BEHR_NO2)
    if DEBUG_LEVEL > 1; fprintf('   Pixel %d of %d\n',a,numel(BEHR_NO2)); end
    IN = inpolygon(lon,lat,corner_lon(:,a),corner_lat(:,a));
    DC8_NO2(a) = nanmean(no2_integrated_col(IN));
end

% Finally, remove pixels that didn't have any aircraft measurements
notnans = ~isnan(DC8_NO2);
omi_lat = omi_lat(notnans); omi_lon = omi_lon(notnans);
BEHR_NO2 = BEHR_NO2(notnans); OMI_NO2 = OMI_NO2(notnans);
DC8_NO2 = DC8_NO2(notnans);
end

