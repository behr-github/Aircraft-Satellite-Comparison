function [ omi_lon_out, omi_lat_out, omi_no2_out, behr_no2_out, air_no2_out, coverage_fraction, quality_flags, db ] = spiral_verification( Merge, Data, timezone, varargin )
%[lon, lat, omi, behr, air, cov_frac] = spiral_verification(Merge,Data,timezone) Compare OMI pixel NO2 values to aircraft spirals.
%
%   This function is based off of the method described in Hains et. al. (J.
%   Geophys. Res. 2010, 115, D05301 doi 10.1029/2009JD012399) to compare
%   OMI NO2 columns to aircraft measurements using full spiral profiles as
%   the means to integrate NO2 concentrations.
%
%   This requires 3 inputs: Merge - a data structure (produced by
%   read_merge_data.m) that contains all the data for one day of flight
%   campaigns. Data is one top-level index of one of the Data data
%   structures output by BEHR_main.m. Timezone must be one of the four
%   standard US timezones (est, cst, mst, or pst) and is needed to compare
%   UTC time to OMI overpass time.
%
%   This outputs pixel longitude and latitude, OMI and BEHR satellite NO2
%   columns for those pixels, the corresponding aircraft inferred columns,
%   and the fraction of spiral measurments that fall within the pixel
%   boundary.
%
%   Parameters:
%
%   profiles: Allows the user to pass either (1) the name of
%   the field containing profile ID numbers in the Merge structure, or(2)
%   an (n x 2) matrix containing the start and stop times (in seconds after
%   midnight UTC) of the periods during the flight when the aircraft is
%   sprialing.  The first is useful if the profile numbers field is not
%   recognized by this function; the second is useful for campaigns such as
%   ARCTAS-CA that do not identify spirals.
%
%   no2field: Defaults to 'NO2_LIF', if this is not the NO2 field, use this
%   parameter to override that.
%
%   radarfield: Defaults to the correct radar altitude field name for a
%   DISCOVER campaign based on the date of the merge file.  Use this field
%   to override that selection.
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
%
%   Josh Laughner <joshlaugh5@gmail.com> 4 Jul 2014

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst',})));
p.addParamValue('profiles',[], @(x) size(x,2)==2 || ischar(x));
p.addParamValue('no2field','NO2_LIF',@isstr);
p.addParamValue('radarfield','',@isstr);
p.addParamValue('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis'})));
p.addParamValue('cloud_frac_max',0.2, @isscalar);
p.addParamValue('rowanomaly','AlwaysByRow',@(x) strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'}));
p.addParamValue('DEBUG_LEVEL',1,@isscalar);

p.parse(Merge,Data,timezone,varargin{:});
pout = p.Results;

% Check that only one element of Data was passed
if numel(Data)>1; error('bdy_layer_verify:DataInput','Only pass one top-level element of the Data structure'); end

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
profiles = pout.profiles;
no2field = pout.no2field;
radarfield = pout.radarfield;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

merge_datenum = datenum(Merge.metadata.date);

if isempty(profiles)
    spiral_mode = 'profnum';
    if merge_datenum >= datenum('07/01/2011') && merge_datenum <= datenum('07/31/2011');
        profnum = Merge.Data.ProfileSequenceNum.Values; % DISCOVER-AQ in Baltimore
    elseif merge_datenum >= datenum('01/16/2013') && merge_datenum <= datenum('02/06/2013');
        profnum = Merge.Data.ProfileNumber.Values; % DISCOVER-AQ in CA
    elseif merge_datenum >= datenum('09/04/2013') && merge_datenum <= datenum('09/29/2013');
        profnum = Merge.Data.ProfileNumber.Values; % DISCOVER-AQ in Texas
    else
        error('sprial_ver:profnum','Profiles not identified. Pass a profile number field name or UTC time ranges as the parameter ''profiles''');
    end
elseif ischar(profiles)
    spiral_mode = 'profnum';
    profnum = eval(sprintf('Merge.Data.%s.Values',profiles));
elseif ismatrix(profiles);
    spiral_mode = 'utcranges';
    Ranges = profiles;
end

if isempty(radarfield)
    if merge_datenum >= datenum('07/01/2011') && merge_datenum <= datenum('07/31/2011');
        radarfield = 'A_RadarAlt';
    elseif merge_datenum >= datenum('01/16/2013') && merge_datenum <= datenum('02/06/2013');
        radarfield = 'Radar_Altitude';
    elseif merge_datenum >= datenum('09/04/2013') && merge_datenum <= datenum('09/29/2013');
        radarfield = 'Radar_Altitude';
    else
        error('sprial_ver:profnum','Radar Alt not identified. Pass a radar altitude field name as the parameter ''radarfield''');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     LOAD DATA     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First handle the aircraft data
[no2, utc, alt, lon, lat] = remove_merge_fills(Merge,no2field);
radar_alt = remove_merge_fills(Merge,radarfield);
pres = remove_merge_fills(Merge,'PRESSURE');
temperature = remove_merge_fills(Merge,'TEMPERATURE');

% Make a composite profile for all data within 3 hours of OMI overpass.
% Find the median of the 10 highest altitude NO2 values (along with their
% associated temperature and pressure values) and append these as the top
% of tropopause value.
tt = utc >= local2utc('10:45',tz) & utc <= local2utc('16:45',tz);
[no2_composite, alt_composite] = bin_rolling_vertical_profile(alt(tt),no2(tt),0.5,0.1);
[alt_composite, no2_composite] = fill_nans(alt_composite, no2_composite); % Remove trailing NaNs and interpolate any internal ones
temp_composite = bin_rolling_vertical_profile(alt(tt), temperature(tt), 0.5, 0.1,'binstart','defined','bincenter',alt_composite);
pres_composite = bin_rolling_vertical_profile(alt(tt), pres(tt), 0.5, 0.1,'binstart','defined','bincenter',alt_composite);

M1 = sortrows([alt_composite', no2_composite', temp_composite', pres_composite']);
top = find(~isnan(M1(:,2)),10,'last');
no2_comp_top_med = median(M1(top,2)); temp_comp_top_med = nanmedian(M1(top,3)); pres_comp_top_med = nanmedian(M1(top,4));
no2_composite(end+1) = no2_comp_top_med; alt_composite(end+1) = 12;
temp_composite(end+1) = temp_comp_top_med; pres_composite(end+1) = pres_comp_top_med;

% Identify all spirals according to the 'profiles' input; reject any
% without a start time between 10:45 and 4:45, that is, about 3 hours on
% either side of the OMI overpass
if strcmp(spiral_mode,'profnum')
    % Get all unique profile numbers and their start times
    unique_profnums = unique(profnum(profnum~=0));
    start_times = zeros(numel(unique_profnums),1);
    for a=1:numel(unique_profnums)
        xx = profnum == unique_profnums(a);
        start_times(a) = min(utc(xx));
    end
    
    % Remove from consideration any profiles with a start time before 10:45
    % am or after 4:45 pm local standard time
    yy = start_times >= local2utc('10:45',tz) & start_times <= local2utc('16:45',tz);
    unique_profnums = unique_profnums(yy); start_times = start_times(yy);
    
    % Save each profile's NO2, altitude, radar altitude, latitude, and
    % longitude as an entry in a cell array
    s = size(unique_profnums);
    no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
    lat_array = cell(s); lon_array = cell(s);
    pres_array = cell(s); temp_array = cell(s);
    for a=1:numel(unique_profnums)
        xx = profnum == unique_profnums(a);
        no2_array{a} = no2(xx);
        alt_array{a} = alt(xx);
        radar_array{a} = radar_alt(xx);
        lat_array{a} = lat(xx);
        lon_array{a} = lon(xx);
        pres_array{a} = pres(xx);
        temp_array{a} = temperature(xx);
    end
elseif strcmp(spiral_mode,'utcranges')
    % Find all the utc start times that are between 10:45 and 4:45 local
    % standard time
    yy = Ranges(:,1) >= local2utc('10:45',tz) & Ranges(:,1) <= local2utc('16:45',tz);
    s = [1,sum(yy)];
    no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
    lat_array = cell(s); lon_array = cell(s);
    pres_array = cell(s); temp_array = cell(s);
    for a=1:s(1)
        xx = utc >= Ranges(a,1) & utc <= Ranges(a,2);
        no2_array{a} = no2(xx);
        alt_array{a} = alt(xx);
        radar_array{a} = radar_alt(xx);
        lat_array{a} = lat(xx);
        lon_array{a} = lon(xx);
        pres_array{a} = pres(xx);
        temp_array{a} = temperature(xx);
    end
end

% Then get pixel information from BEHR: column NO2, pixel corners, etc. We
% will not consider any pixels that were affected by the row anomaly, have
% too great a cloud fraction, etc.

Data.Areaweight = ones(size(Data.Longitude));
if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
Data2 = omi_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
xx = Data2.Areaweight > 0;

omi_lat = Data2.Latitude(xx); omi_lon = Data2.Longitude(xx);
corner_lat = Data2.Latcorn(:,xx); corner_lon = Data2.Loncorn(:,xx);
behr_no2 = Data2.BEHRColumnAmountNO2Trop(xx);
omi_no2 = Data2.ColumnAmountNO2Trop(xx);
TropopausePres = Data2.TropopausePressure(xx);
vza = Data2.ViewingZenithAngle(xx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     MATCH PIXELS AND SPIRALS     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrices to hold the values for each pixel measured.  We will
% over estimate the number of pixels that can be compared, assuming a
% maximum of 4 per spiral.  This prevents MATLAB from needing to rewrite
% the variable to different memory as it grows.
n = numel(no2_array); P = 0; % P will keep track of which index to put the next entry into
omi_lon_out = -9e9*ones(n*4,1);
omi_lat_out = -9e9*ones(n*4,1);
omi_no2_out = -9e9*ones(n*4,1);
behr_no2_out = -9e9*ones(n*4,1);
air_no2_out = -9e9*ones(n*4,1);
coverage_fraction = -9e9*ones(n*4,1);
quality_flags = uint16(zeros(n*4,1));

db.loncorn = -9e9*ones(4,n*4);
db.latcorn = -9e9*ones(4,n*4);

% Also, define Avogadro's number and gas constant;
Av = 6.022e23; % molec. / mol
R = 8.314e4; % (hPa * cm^3) / (mol * K)

% Now iterate through each profile, and find all pixels which that profile
% overlaps. As Hains did, we will require that the pixel contains
% measurements between 0-3 km, specifically (in addition to what was
% described in her paper), we will require that there be 20 measurements
% below 3 km for that pixel, and that the pixel's VZA be less than 60
% degrees.
for p=1:n
    if min(radar_array{p}) > 0.5
        continue % As per Hains et. al., the bottom of the profile must reach down to at least 500 m above the surface
    else
        % Bin the NO2 data by altitude.  The center of the bottom bin will
        % be at the lowest altitude that has NO2 data
        [no2bins, altbins] = bin_rolling_vertical_profile(alt_array{p}, no2_array{p}, 0.5, 0.1,'binstart','center');
        [altbins, no2bins] = fill_nans(altbins,no2bins);
        [tempbins, temp_altbins] = bin_rolling_vertical_profile(alt_array{p}, temp_array{p}, 0.5, 0.1, 'binstart', 'defined', 'bincenters', altbins);
        [presbins, pres_altbins] = bin_rolling_vertical_profile(alt_array{p}, pres_array{p}, 0.5, 0.1, 'binstart', 'defined', 'bincenters', altbins);
        % Figure out how far above the surface the plane was at its lowest,
        % then add an extra bin at the bottom that has an NO2 concentration
        % equal to the median of the lowest 10 NO2 measurements. Also, add
        % an extra bin at the top that will hold the value extrapolated to
        % the tropopause pressure.
        M = sortrows([alt_array{p}', no2_array{p}', radar_array{p}', pres_array{p}', temp_array{p}']);
        xx = find(~isnan(M(:,2)),10,'first'); zz = find(~isnan(M(:,2)),10,'last');
        med_no2 = median(M(xx,2)); top_med_no2 = median(M(zz,2));
        med_pres = nanmedian(M(xx,4)); top_med_pres = nanmedian(M(zz,4));
        med_temp = nanmedian(M(xx,5)); top_med_temp = nanmedian(M(zz,5));
        med_radar_alt = nanmedian(M(xx,3)); surface_alt = altbins(1)-med_radar_alt;
        
        % Initialize the qualtity flag for this pixel
        q_flag = uint16(0);
        
        % If the top median no2 value is < 100 pptv, then it's safe to
        % assume that the profile has crossed the boundary layer and is
        % sampling the free troposphere.  If not, then extrapolating that
        % value to the tropopause will grossly overestimate the total
        % column.  In the latter case, append the composite profile on top
        % of the current one.
        if top_med_no2 < 100;
            no2bins = [med_no2, no2bins, top_med_no2]; altbins = [surface_alt, altbins, 12];
            tempbins = [med_temp, tempbins, top_med_temp];
            presbins = [med_pres, presbins, top_med_pres];
        else
            % Get the altitude of the highest bin in the profile and find
            % all bins in the composite profile above that
            xx = alt_composite > max(altbins);
            no2bins = [med_no2, no2bins, no2_composite(xx)];
            altbins = [surface_alt, altbins, alt_composite(xx)];
            tempbins = [med_temp, tempbins, temp_composite(xx)];
            presbins = [med_pres, presbins, pres_composite(xx)];
            
            % Set the 3rd bit of the quality flag to 1 to indicate that a
            % composite column was appended
            q_flag = bitset(q_flag,3,1);
        end
        
        % Remove pixels not overlapping this profile
        if DEBUG_LEVEL > 0; fprintf('   Removing pixels outside profile track\n'); end
        lon_logical = corner_lon > min(lon_array{p}) & corner_lon < max(lon_array{p});
        lat_logical = corner_lat > min(lat_array{p}) & corner_lat < max(lat_array{p});
        latlon_logical = any(lon_logical) & any(lat_logical);
        
        omi_lat_p = omi_lat(latlon_logical); omi_lon_p = omi_lon(latlon_logical);
        corner_lat_p = corner_lat(:,latlon_logical); corner_lon_p = corner_lon(:,latlon_logical);
        behr_no2_p = behr_no2(latlon_logical); omi_no2_p = omi_no2(latlon_logical);
        TropopausePres_p = TropopausePres(latlon_logical); vza_p = vza(latlon_logical);
        
        for o=1:numel(behr_no2_p);
            no2_3km = no2_array{p}(alt_array{p}<3);
            lon_3km = lon_array{p}(alt_array{p}<3);
            lat_3km = lat_array{p}(alt_array{p}<3);
            IN_3km = inpolygon(lon_3km, lat_3km, corner_lon_p(:,o), corner_lat_p(:,o));
            if sum(~isnan(no2_3km(IN_3km)))<20
                continue % Pixel must have 20 valid measurments between 0-3 km altitude (good sampling of boundary layer)
            elseif vza_p(o) > 60
                continue % Pixel viewing zenith angle must not be greater than 60 degress - pixels on the edge of the track are very wide and the spiral is less likely to be representative.
            else
                P = P+1;
                % Insert the OMI tropopause altitude as the final altitude
                % bin, then interpolate the NO2, temperature, and pressure
                % values and integrate
                TP_alt = -log(TropopausePres_p(o)/1013.25)*7.4;
                altbins(end) = TP_alt;
                
                % Interpolate the NO2, temperature, and pressure data
                dz = 1; % integration segments in meters
                alt_profile = altbins(1):(dz/1000):altbins(end);
                no2_profile = interp1(altbins,no2bins,alt_profile,'linear');
                temp_profile = interp1(altbins,tempbins,alt_profile,'linear');
                pres_profile = exp(interp1(altbins,log(presbins),alt_profile,'linear')); % Linearly interpolate ln(P) since that is what depends linearly on altitude
                
                % Carry out the numerical integration
                no2_column = 0;
                for z=1:numel(alt_profile)
                    P_z = pres_profile(z); T = temp_profile(z); no2_z = no2_profile(z);
                    conc_NO2 = (Av * P_z * no2_z * 1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
                    no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
                end
                
                % If any bits in the quality flag are set, set the summary
                % bit
                if any(q_flag); q_flag = bitset(q_flag,1,1); end
                
                % Save the results for this pixel
                omi_lon_out(P) = omi_lon_p(o);
                omi_lat_out(P) = omi_lat_p(o);
                omi_no2_out(P) = omi_no2_p(o);
                behr_no2_out(P) = behr_no2_p(o);
                air_no2_out(P) = no2_column;
                quality_flags(P) = q_flag;
                
                db.latcorn(:,P) = corner_lat_p(:,o);
                db.loncorn(:,P) = corner_lon_p(:,o);
                
                IN_all = inpolygon(lon_array{p}, lat_array{p}, corner_lon_p(:,o), corner_lat_p(:,o));
                coverage_fraction(P) = sum(IN_all)/numel(no2_array{p});
            end
        end
    end
end
% Clean up the output variables
fills = air_no2_out == -9e9;
omi_lon_out = omi_lon_out(~fills);
omi_lat_out = omi_lat_out(~fills);
omi_no2_out = omi_no2_out(~fills);
behr_no2_out = behr_no2_out(~fills);
air_no2_out = air_no2_out(~fills);
coverage_fraction = coverage_fraction(~fills);
quality_flags = quality_flags(~fills);

db.latcorn = db.latcorn(:,~fills);
db.loncorn = db.loncorn(:,~fills);
end

