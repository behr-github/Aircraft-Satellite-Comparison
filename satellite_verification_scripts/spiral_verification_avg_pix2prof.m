function [ prof_lon_out, prof_lat_out, omi_no2_out, behr_no2_out, air_no2_out, db ] = spiral_verification_avg_pix2prof( Merge, Data, timezone, varargin )
%[lon, lat, omi, behr, air, cov_frac] = spiral_verification(Merge,Data,timezone) Compare OMI pixel NO2 values to aircraft spirals.
%
%   This function is based off of the method described in Hains et. al. (J.
%   Geophys. Res. 2010, 115, D05301 doi 10.1029/2009JD012399) and Bucesela
%   et. al. (J. Geophys. Res. 2008, 113, D16S31 doi:10.1029/2007/D008838)
%   to compare OMI NO2 columns to aircraft measurements using full spiral
%   profiles as the means to integrate NO2 concentrations.  The profile
%   data is binned by pressure bins as in the Bucsela article.
%   Extrapolation above and below the spiral is usually handled as in Hains
%   et. al.; the median lowest/highest 10 NO2 measurements are extrapolated
%   as a constant.  As in Hains, the surface altitude is determined as the
%   difference between the median of the lowest 10 GPS altitudes and radar
%   altitudes. (Pressure altitude was not used because of potential
%   variability with local conditions).  Unlike Hains, the tropopause
%   pressure was not fixed, but taken from the OMI pixel being validated
%   against.
%
%   Unlike spiral_verification which averages all profiles to the pixel
%   they cross into, this function averages all measurements for the pixels
%   that cover a particular profile together.
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
%   The sixth output is a structure containing a wide range of ancillary
%   data that might be valuable for debugging, such as individual pixel
%   measurements, cloud fractions, etc. All fields contain cell arrays
%   where each cell corresponds to a profile. Most fields are self
%   explanatory, but two require some additional explanation.
%       The quality_flags field is a 16-bit integer that is a quality flag
%   for the column, similar to the vcdQualityFlag and XTrackFlag in the OMNO2
%   product.  Use bitget() to check the status of individual bits; the
%   meaning of each bit is specified here:
%       1st: Summary, set to 1 if any flags are set.
%       2nd: Reserved as a second summary bit against future need.
%       3rd: Indicates the column top was derived from a daily composite
%           rather than extrapolation of the top median NO2 values.
%       4th: Indicates that < 10% of the data points in the profile had NO2
%           data
%       5th: Indicates that the composite profile had < 10% of the data
%           points valid
%       6th: Indicates that radar values from higher in the column were
%           used to calculate the surface pressure.
%       7th: Indicates that the GLOBE terrain database was used to find the
%           surface pressure because no radar data was available
%       8th: No data points fall within the time frame of +/- 3 hr from OMI
%           overpass
%       9th: No BEHR data for this swath (probably used an OMI_SP only
%           file)
%       10th: A warning that the flight path is considered to be in
%           multiple time zones.
%       11-15: Unused
%       16th: Set if the column was skipped due to < 1% valid
%           NO2, pressure, or temperature data
%
%       The reject field contains more specific information on why a
% profile or pixel was rejected.  The difference between this and the
% quality_flags field is that the quality_flags field tries to represent
% problems with the underlying data (and generally involves the whole merge
% dataset), while the reject field describes reasons that a profile or the
% pixels impinging on it were not considered for comparison.  Generally, if
% a cell contains a matrix, each entry represents one pixel impinging on
% the profile. It is an 8-bit integer; the bits have the following
% meanings:
%       1st: Profile did not reach minRadarAlt of the ground (from Hains
%           et. al., recommended 0.5 km) 
%       2nd: Entire profile was NaNs 
%       3rd: Less than numBLpoints data points in the lowest 3 km (from Hains et.
%          al., recommended number is 20).
%       4th: No valid pixels (clouds, row anomaly) overlap profile
%       5th: VZA for this pixel > 60 deg.
%       6th: Profile height was too little (max(alt) - min(alt) <
%          min_height). The required height is set with the parameter
%          'min_height'
%
%   Parameters:
%
%   campaign_name: a string identifying the campaign, used to retrieve the
%   appropriate field names. See merge_field_names.m in the Utils/Constants
%   folder.  The string need only contain an identifiable campaign name,
%   note that the Discover campaigns need the string 'discover' plus the
%   state abbreviation.  An empty string will bypass this, note that
%   no2field, altfield, radarfield, and aerfield will need to be set then
%   or an error will be thrown.
%
%   profiles: Allows the user to pass either (1) the name of the field
%   containing profile ID numbers in the Merge structure, (2) an (n x 2)
%   matrix containing the start and stop times (in seconds after midnight
%   UTC) of the periods during the flight when the aircraft is sprialing,
%   or The first is useful if the profile numbers field is not recognized
%   by this function; the second is useful for campaigns such as ARCTAS-CA
%   that do not identify spirals.
%
%   profnums: A list of specific profile numbers to examine, usually used
%   to examine only certain profiles that have been pre-selected, for e.g.
%   their aerosol characteristics.  Defaults to an empty matrix, which
%   means all profiles are included.
%
%   behrfield: What field to use for BEHR NO2 data. Defaults to
%   BEHRColumnAmountNO2Trop, but can be reset to, for example, use
%   reprocessed columns using in-situ AMFs.
%
%   starttime: Profiles must have a start time later that this. Pass as a
%   string using military time, e.g. 16:00 instead of 4:00 pm.  This is
%   always in local standard time.  Set to 10:45 by default.  Allows user
%   to restrict aircraft data to times near satellite overpass.
%
%   endtime: Profiles must have a start time before this. See starttime for
%   more details.  Set to 16:45 by default.
%
%   no2field: If unspecified, passed an empty string, or set to 'lif' will
%   default to the LIF (our data) field, field name determined by calling
%   merge_field_names with the campaign name specified. If set to 'cl',
%   will use the chemiluminescence data - NCAR if available, NOAA
%   otherwise.. Otherwise will use the field specified by the string given.
%
%   altfield: If unspecified or set to 'gps' uses the GPS altitude field
%   for the campaign, if set to 'pressure', uses the pressure altitude
%   field. Otherwise will use the field specified by the string given.
%
%   radarfield: If unspecified uses the default radar altitude field for
%   the given campaign, otherwise uses the field specified by this input.
%
%   presfield: Defaults to PRESSURE.
%
%   tempfield: Defaults to TEMPERATURE.
%
%   aerfield: If set to an empty string (default) will try to figure out
%   the aerosol extinction field based on the dates given. If set to 0,
%   will not return aerosol data. If set to a string, will use that string
%   as the field name.
%
%   cloud_product: Which cloud product (omi or modis or rad) to use in
%   rejecting pixels.  Defaults to omi.
%
%   cloud_frac_max: The maximum allowed geometric cloud fraction.  Defaults
%   to 0.2; recommended value for use with MODIS cloud product is 0.
%
%   rowanomaly: The method of rejecting pixels based on row anomaly.
%   Defaults to 'AlwaysByRow'.  See omi_pixel_reject or omi_rowanomaly for
%   more information on the possible choices ('AlwaysByRow', 'RowsByTime',
%   'XTrackFlags', and 'XTrackFlagsLight').
%
%   min_height: The minimum difference between the lowest and highest
%   points in the profile. Defaults to 0, i.e. any profile height.
%
%   numBLpoints: The minimum number of data points (valid, not NaNs) in the
%   lowest 3 km of the atmosphere.  Hains et. al. recommends 20, which is
%   set as the default.
%
%   minRadarAlt: The altitude (in km) above the ground which the plane must
%   go below for the profile to be used. Hains et. al. recommend 0.5 km
%   which is the default.
%
%   DEBUG_LEVEL: The level of output messages to write; 0 = none, 1 =
%   normal; 2 = verbose, 3 = (reserved), 4 = plot NO2 profiles colored by
%   source - intrinsic, extrapolation, composite
%
%   clean: Setting this to 0 will not remove any pixels with a fill value -
%   useful only for debugging why a pixel is rejected.
%
%   Josh Laughner <joshlaugh5@gmail.com> 25 Jul 2014


% TODO: add a filter for profiles that can be controlled by optional
% inputs, e.g. how close to the surface, how much breadth of altitude
% covered, max altitude, lat/lon spread.
p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst','auto'})));
p.addParameter('behrfield','BEHRColumnAmountNO2Trop',@isstr);
p.addParameter('starttime','10:45',@isstr);
p.addParameter('endtime','16:45',@isstr);
p.addParameter('campaign_name','',@isstr);
p.addParameter('profiles',[], @(x) size(x,2)==2 || ischar(x));
p.addParameter('profnums',[], @(x) (isvector(x) || isempty(x) || size(x,2)==2));
p.addParameter('no2field','',@isstr);
p.addParameter('conv_fact',1e-12,@isscalar);
p.addParameter('altfield','',@isstr);
p.addParameter('radarfield','',@isstr);
p.addParameter('presfield','PRESSURE',@isstr);
p.addParameter('tempfield','TEMPERATURE',@isstr)
p.addParameter('aerfield','', @(x) (ischar(x) || x==1 || x==0) );
p.addParameter('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis','rad'})));
p.addParameter('cloud_frac_max',0.2, @isscalar);
p.addParameter('rowanomaly','AlwaysByRow',@(x) any(strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'})));
p.addParameter('min_height',0,@isscalar);
p.addParameter('numBLpoints',20,@isscalar);
p.addParameter('minRadarAlt',0.5,@isscalar);
p.addParameter('DEBUG_LEVEL',1,@isscalar);
p.addParameter('clean',1,@isscalar);

p.parse(Merge,Data,timezone,varargin{:});
pout = p.Results;

% Check that only one element of Data was passed
if numel(Data)>1; error('bdy_layer_verify:DataInput','Only pass one top-level element of the Data structure'); end

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
behrfield = pout.behrfield;
starttime = pout.starttime;
endtime = pout.endtime;
profiles = pout.profiles;
user_profnums = pout.profnums;
campaign_name = pout.campaign_name;
no2field = pout.no2field;
conv_fact = pout.conv_fact;
altfield = pout.altfield;
radarfield = pout.radarfield;
presfield = pout.presfield;
aerfield = pout.aerfield;
Tfield = pout.tempfield;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;
min_height = pout.min_height;
numBLpoints = pout.numBLpoints;
minRadarAlt = pout.minRadarAlt;
DEBUG_LEVEL = pout.DEBUG_LEVEL;
clean_bool = pout.clean;

merge_datenum = datenum(Merge.metadata.date);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Many of the errors in this script have not been updated yet to use the
% new error class.  
E = JLLErrors;


% As long as the campaign name isn't blank, try to get the appropriate
% field names for the given campaign. Otherwise check that the four merge
% fields we need are set and if not, error out.
if ~isempty(campaign_name);
    FieldNames = merge_field_names(campaign_name);
elseif isempty(no2field) || isempty(altfield) || isempty(aerfield) || isempty(radarfield)
    error(E.badinput('If no campaign is to be specified (using the parameter ''campaign_name'') then parameters no2field, altfield, radarfield, and aerfield must not be empty strings.'));
end

% Deal with the profiles or UTC range input. Set both the spiral_mode
% string (which determines whether to use profile numbers or UTC ranges
% later on in this function) and reads in the profile numbers or UTC ranges.
if isempty(profiles)
    spiral_mode = 'profnum';
    if isempty(FieldNames.profile_numbers)
        error(E.callError('profile_mode','The profile numbers field is not defined for this campaign; be sure to set the calling function to use UTC ranges'));
    end
    profnum = Merge.Data.(FieldNames.profile_numbers).Values;
elseif ischar(profiles)
    spiral_mode = 'profnum';
    profnum = Merge.Data.(profiles).Values;
elseif ismatrix(profiles);
    spiral_mode = 'utcranges';
    Ranges = profiles;
end

% Now for each field name not given by the user, set it from the
% merge field names retrieved.
if isempty(no2field) || strcmpi(no2field,'lif')
    no2field = FieldNames.no2_lif;
elseif strcmpi(no2field,'cl')
    no2field = FieldNames.no2_ncar;
end

if isempty(altfield) || strcmpi(altfield,'gps')
    altfield = FieldNames.gps_alt;
elseif strcmpi(altfield,'pressure')
    altfield = FieldNames.pressure_alt;
end

if isempty(radarfield)
    radarfield = FieldNames.radar_alt;
end

if isempty(aerfield) % This will not override a 0 passed as this field.
    aerfield = FieldNames.aerosol_extinction;
end

% Finally check that user_profnums fullfills the requirements for the
% profile mode, be it profile numbers or UTC ranges.
if ~isempty(user_profnums)
    if strcmp(spiral_mode,'profnum')
        if ~isvector(user_profnums)
            error(E.badinput('User defined profile numbers must be a vector.'));
        end
    elseif strcmp(spiral_mode,'utcranges')
        if size(user_profnums,2) ~= 2
            error(E.badinput('User defined UTC ranges must be in an n x 2 matrix'));
        end
    end
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     LOAD DATA     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First handle the aircraft data
[no2, utc, alt, lon, lat] = remove_merge_fills(Merge,no2field,'alt',altfield);
no2(no2<0) = NaN; % Must remove any negative values from consideration because they will return imaginary components during the log-log interpolation
lon(isnan(lat))=NaN; lat(isnan(lon))=NaN; % Make any points that are NaNs in one fields also so in the other
radar_alt = remove_merge_fills(Merge,radarfield,'alt',altfield);
pres = remove_merge_fills(Merge,presfield,'alt',altfield);
temperature = remove_merge_fills(Merge,Tfield,'alt',altfield);
altfill = eval(sprintf('Merge.Data.%s.Fill',altfield));
alt(alt==altfill) = NaN; % Switching to GPS altitude gave fill values for altitude.  These must be removed.
if aerfield ~= 0
    aer_data = remove_merge_fills(Merge,aerfield);
end

% The variable that holds the quality flags
q_base = uint16(0);

% If the timezone was set to "auto," calculate the difference from UTC
% based on the mean longitude. This will produce a vector of the same
% length as the utc and lon variables
if strcmpi(tz,'auto')
    tz = round(lon/15);
    
    % If all elements in the timezone vector aren't the same, set the 10th
    % bit in the quality flag
    if ~all(tz==tz(1))
        q_base = bitset(q_base,10,1);
    end
end

% Now, use the timezone (entered or calculated) to produce a new vector of
% times that reflect the local time at each point
local_times = utc2local_sec(utc,tz);

% Check what percentage of values were fill values, if >99% are fills for
% data, temperature, or pressure, return NaNs and set the largest bit on
% the quality flag to 1 (as well as the summary flag). If <99% but >90%,
% set the 5th bit to 1 as a warning.
percent_nans = sum(isnan(no2))/numel(no2);
percent_nans_P = sum(isnan(pres))/numel(pres);
percent_nans_T = sum(isnan(temperature))/numel(temperature);
if percent_nans > 0.99 || percent_nans_P > 0.99 || percent_nans_T > 0.99;
    if DEBUG_LEVEL > 1 && percent_nans > 0.99;
        fprintf('%s had < 1%% of NO2 values valid\n',datestr(merge_datenum,2));
    end
    if DEBUG_LEVEL > 1 && percent_nans_P > 0.99;
        fprintf('%s had < 1%% of pressure values valid\n',datestr(merge_datenum,2));
    end
    if DEBUG_LEVEL > 1 && percent_nans_T > 0.99;
        fprintf('%s had < 1%% of temperature values valid\n',datestr(merge_datenum,2));
    end
    % We must skip this merge file if there is no NO2, pressure, or
    % temperature data
    prof_lon_out = NaN;
    prof_lat_out = NaN;
    omi_no2_out = NaN;
    behr_no2_out = NaN;
    air_no2_out = NaN;
    
    db.all_omi = {NaN};
    db.all_behr = {NaN};
    db.quality_flags = {uint16(2^15+1)};
    db.coverage_fraction = {NaN};
    db.dist_vectors = {NaN};
    db.loncorn = {NaN(4,1)};
    db.latcorn = {NaN(4,1)};
    db.strat_NO2 = {NaN};
    db.modis_cloud = {NaN};
    if strcmp(spiral_mode,'profnum'); db.profnums = {NaN};
    else db.profnums = {NaN(1,2)};
    end
    db.reject = {uint8(2)};
    db.lon_3km = {NaN};
    db.lat_3km = {NaN};
    if aerfield ~= 0; 
        db.aer_max_out = {NaN}; 
        db.aer_int_out = {NaN};
        db.aer_quality = {NaN};
    end
    db.column_error = {NaN};
else
    if percent_nans > 0.9;
        warning('Merge file for %s has %.0f%% NO2 values as NaNs',datestr(merge_datenum,2),percent_nans*100);
        q_base = uint16(bitset(q_base,5,1));
    end
    
    
    % Make a composite profile for all data within 3 hours of OMI overpass
    % - this will look for points based on their local time.
    % Find the median of the 10 highest altitude NO2 values (along with their
    % associated temperature and pressure values) and append these as the top
    % of tropopause value.
    tt = local_times >= local2utc('10:45',0) & local_times <= local2utc('16:45',0);
    if sum(tt) == 0 % If no points fall within the time frame, return NaNs and exit
        prof_lon_out = NaN;
        prof_lat_out = NaN;
        omi_no2_out = NaN;
        behr_no2_out = NaN;
        air_no2_out = NaN;
        
        db.all_omi = {NaN};
        db.all_behr = {NaN};
        q_base = bitset(q_base,1,1); db.quality_flags = {bitset(q_base,8,1)};
        db.coverage_fraction = {NaN};
        db.dist_vectors = {NaN};
        db.loncorn = {NaN(4,1)};
        db.latcorn = {NaN(4,1)};
        db.strat_NO2 = {NaN};
        db.modis_cloud = {NaN};
        if strcmp(spiral_mode,'profnum'); db.profnums = {NaN};
        else db.profnums = {NaN(1,2)};
        end
        db.reject = {uint8(2)};
        db.lon_3km = {NaN};
        db.lat_3km = {NaN};
        if aerfield ~= 0; 
            db.aer_max_out = {NaN}; 
            db.aer_int_out = {NaN};
            db.aer_quality = {NaN};
        end
        db.column_error = {NaN};
        return
    end
    [no2_composite, pres_composite] = bin_omisp_pressure(pres(tt),no2(tt));
    % Calculate the uncertainty in the composite bins as the std. error
    [~,~, no2_composite_stderr] = bin_omisp_pressure(pres(tt),no2(tt),'mean');
    temp_composite = bin_omisp_pressure(pres(tt),temperature(tt));
    
    % The top bin of bin_omisp_pressure (200 hPa) is right around the normal
    % boundary of the troposphere, 12 km.  If no NO2 data is available, (i.e.
    % that bin has a value of NaN), then we'll extrapolate the median of the
    % top ten NO2 measurements to that bin.  In the loop itself, we'll adjust
    % for the changing tropopause pressure.
    if isnan(no2_composite(end))
        M1 = sortrows([pres', no2', temperature']);
        top = find(~isnan(M1(:,2)),10,'first'); % Since lower pressure = higher altitude, we want the first 10 NO2 measurements when sorted by pressure.
        no2_comp_top_med = median(M1(top,2)); temp_comp_top_med = nanmedian(M1(top,3)); pres_comp_top_med = nanmedian(M1(top,1));
        
        % Assume the error is the standard error of the top ten points
        no2_comp_top_med_stderr = std(M1(top,2))/sqrt(10);
        if no2_comp_top_med < 3; 
            no2_comp_top_med = 1.5; 
            if DEBUG_LEVEL > 0; fprintf('Composite profile top < LoD\n'); end
        end
        no2_composite(end) = no2_comp_top_med;
        no2_composite_stderr(end) = no2_comp_top_med_stderr;
        
        % Temperature should have a linear relationship to altitude
        % (i.e., log(pressure)) up to the tropopause; therefore extrapolate
        % based on the fit.
        nans = isnan(temp_composite);
        temp_composite(nans) = interp1(log(pres_composite(~nans)), temp_composite(~nans), log(pres_composite(nans)),'linear','extrap');
        
        % Do any interpolation in log-log space (per. Bucsela et. al. J.
        % Geophys. Res. 2008, 113, D16S31) since pressure has an exponential
        % relation to altitude.  Log-log space would assume that NO2 also has
        % an exponential relationship to altitude.
        [~, no2_tmp, no2_tmperr] = fill_nans(log(pres_composite),log(no2_composite),log(no2_composite_stderr),'noclip');
        no2_composite = exp(no2_tmp);
        no2_composite_stderr = exp(no2_tmperr);    
        
    end
    
    % Identify all spirals according to the 'profiles' input; reject any
    % without a start time between 10:45 and 4:45, that is, about 3 hours on
    % either side of the OMI overpass
    if strcmp(spiral_mode,'profnum')
        % Get all unique profile numbers and their start times
        unique_profnums = unique(profnum(profnum~=0));
        start_times = zeros(numel(unique_profnums),1);
        for a=1:numel(unique_profnums)
            % If we're using a vector of timezones, get the most common
            % timezone from the profile - we'll treat the whole profile as
            % belonging to that timezone.  If timezone has been set
            % manually, then just use that time zone to convert the start
            % time of the profile. (obviously) Save the local start time.
            xx = profnum == unique_profnums(a);
            if ismatrix(tz) && isnumeric(tz)
                % Case where we're using a vector of timezones
                mct = mode(tz(xx));
                start_times(a) = utc2local_sec(min(utc(xx)),mct);
            elseif ischar(tz)
                start_times(a) = utc2local_sec(min(utc(xx)),tz);
            else
                error(E.callError('tz_not_recognized','Cannot understand the format of timezones'));
            end
        end
        
        % Remove from consideration any profiles with a start time before 10:45
        % am or after 4:45 pm local standard time
        yy = start_times >= local2utc(starttime,0) & start_times <= local2utc(endtime,0);
        unique_profnums = unique_profnums(yy); start_times = start_times(yy);
        
        % If the user passed a list of profile numbers, remove any profile
        % numbers that do not match the list provided.
        if ~isempty(user_profnums)
            tmp = true(size(unique_profnums));
            for a = 1:numel(unique_profnums)
                tmp(a) = any(unique_profnums(a)==user_profnums);
            end
            unique_profnums(~tmp) = [];
        end
        
        % Save each profile's NO2, altitude, radar altitude, latitude, and
        % longitude as an entry in a cell array
        s = size(unique_profnums);
        no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
        lat_array = cell(s); lon_array = cell(s); profnum_array = cell(s);
        pres_array = cell(s); temp_array = cell(s); aer_array = cell(s);
        for a=1:numel(unique_profnums)
            xx = profnum == unique_profnums(a);
            no2_array{a} = no2(xx);
            alt_array{a} = alt(xx);
            radar_array{a} = radar_alt(xx);
            lat_array{a} = lat(xx);
            lon_array{a} = lon(xx);
            pres_array{a} = pres(xx);
            temp_array{a} = temperature(xx);
            if aerfield ~= 0; aer_array{a} = aer_data(xx); end
            profnum_array{a} = unique_profnums(a);
        end
    elseif strcmp(spiral_mode,'utcranges')
        % Find all the utc start times that are between within the
        % specified range of local times.  Go through each range, find the
        % data points that correspond to it, get the most common timezone,
        % use that to set whether to include that range or not. Also, check
        % the "user_profnums" variable which will have specific UTC ranges
        % to allow
        yy = false(size(Ranges,1),1);
        for a=1:size(Ranges,1)
            if ~isempty(user_profnums) && ~(containedin(Ranges(a,1),user_profnums(:,1)) && containedin(Ranges(a,2),user_profnums(:,2)))
                % If the UTC range is not one defined in the user_profnums
                % (which is really user_ranges in this case, but it was
                % called profnums first), skip the rest of this loop, which
                % will not set yy(a) to true thus skipping this range.
                continue
            end
            tz_ind = utc >= Ranges(a,1) & utc <= Ranges(a,2);
            mct = mode(tz(tz_ind));
            range_start_local = utc2local_sec(Ranges(a,1),mct);
            yy(a) = range_start_local >= local2utc(starttime,0) && range_start_local <= local2utc(endtime,0);
        end
        %yy = Ranges(:,1) >= local2utc(starttime,tz) & Ranges(:,1) <= local2utc(endtime,tz);
        ranges_in_time = Ranges(yy,:);
        s = [1,sum(yy)];
        no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
        lat_array = cell(s); lon_array = cell(s); aer_array = cell(s);
        pres_array = cell(s); temp_array = cell(s); profnum_array = cell(s);
        for a=1:s(2)
            xx = utc >= ranges_in_time(a,1) & utc <= ranges_in_time(a,2);
            no2_array{a} = no2(xx);
            alt_array{a} = alt(xx);
            radar_array{a} = radar_alt(xx);
            lat_array{a} = lat(xx);
            lon_array{a} = lon(xx);
            pres_array{a} = pres(xx);
            temp_array{a} = temperature(xx);
            if aerfield ~= 0; aer_array{a} = aer_data(xx); end
            profnum_array{a} = ranges_in_time(a,:);
        end
    end
    
    % Then get pixel information from BEHR: column NO2, pixel corners, etc. We
    % will not consider any pixels that were affected by the row anomaly, have
    % too great a cloud fraction, etc.
    
    Data.Areaweight = ones(size(Data.Longitude));
    if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
    try % Handles the case of both BEHR files and SP-only files
        Data2 = omi_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
    catch err;
        if strcmp(err.identifier,'MATLAB:nonExistentField');
            Data2 = omi_sp_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
        else
            rethrow(err)
        end
    end
    xx = Data2.Areaweight > 0;
    
    omi_lat = Data2.Latitude(xx); omi_lon = Data2.Longitude(xx);
    corner_lat = Data2.Latcorn(:,xx); corner_lon = Data2.Loncorn(:,xx);
    omi_no2 = Data2.ColumnAmountNO2Trop(xx);
    TropopausePres = Data2.TropopausePressure(xx);
    vza = Data2.ViewingZenithAngle(xx);
    TerrainPres = Data2.TerrainPressure(xx); %load terrain pressure (in hPa)
    % Try to load the BEHR column, if it is a field in the Data structure
    try
        behr_no2 = Data2.(behrfield)(xx);
    catch err
        % If not, fill the imported variable with fill values
        if strcmp(err.identifier,'MATLAB:nonExistentField')
            if DEBUG_LEVEL > 0; fprintf('    No BEHR data for this swath\n'); end
            q_base = bitset(q_base,9,1);
            behr_no2 = -127*ones(size(xx));
        else
            rethrow(err)
        end
    end
    % Try to load the MODIS cloud fraction - not needed within this script,
    % but included in the output structure "db" to examine if cloud
    % fraction has an impact on the retrieval
    try
        modis_cloud = Data2.MODISCloud(xx);
    catch err
        % If "MODISCloud" is not a field, or if it only contains a single
        % value of "0" or is empty (thus the field was created but never
        % had anything entered), fill with import variable with fill
        % values.
        if strcmp(err.identifier,'MATLAB:nonExistentField') || (strcmp(err.identifier,'MATLAB:badsubscript') && isempty(Data2.MODISCloud)) || (strcmp(err.identifier,'MATLAB:badsubscript') && numel(Data2.MODISCloud) == 1 && Data2.MODISCloud == 0)
            if DEBUG_LEVEL > 0; fprintf('    No MODIS data for this swath\n'); end
            modis_cloud = -127*ones(size(xx));
        else
            rethrow(err)
        end
    end
    % Extra fields carried through for curiosity; this is used to calculate
    % stratospheric NO2
    total_omi_no2 = Data2.ColumnAmountNO2(xx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     MATCH PIXELS AND SPIRALS     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize matrices to hold the values for each pixel measured.  
    n = numel(no2_array);
    
    prof_lon_out = -9e9*ones(n,1);
    prof_lat_out = -9e9*ones(n,1);
    omi_no2_out = -9e9*ones(n,1);
    behr_no2_out = -9e9*ones(n,1);
    air_no2_out = -9e9*ones(n,1);
    
    
    
    db.all_omi = cell(n,1);
    db.all_behr = cell(n,1);
    db.quality_flags = cell(n,1);
    db.coverage_fraction = cell(n,1);
    db.dist_vectors = cell(n,1);
    db.loncorn = cell(n,1);
    db.latcorn = cell(n,1);
    db.strat_NO2 = cell(n,1);
    db.modis_cloud = cell(n,1);
    db.profnums = cell(n,1);
    db.reject = cell(n,1);
    db.lon_3km = cell(n,1);
    db.lat_3km = cell(n,1);
    if aerfield ~= 0; 
        db.aer_max_out = cell(n,1); 
        db.aer_int_out = cell(n,1);
        db.aer_quality = cell(n,1);
    end
    db.column_error = cell(n,1);
    
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
        % Initialize the qualtity flag for this profile
        q_flag = q_base;
        
        prof_reject = uint8(0);
        
        lontest = lon_array{p}; lontest(isnan(lontest)) = [];
        lattest = lat_array{p}; lattest(isnan(lattest)) = [];
        % If the profile does not go below 500 m, skip it anyway (per
        % Hains, require for good BL sampling)
        if min(radar_array{p}) > minRadarAlt
            prof_reject = bitset(prof_reject,1); 
            db.reject{p} = prof_reject;
            continue
            % If the entire profile was fill values (i.e. instrument trouble) obviously we have to skip this profile.
        elseif all(isnan(no2_array{p}));
            prof_reject = bitset(prof_reject,2); 
            db.reject{p} = prof_reject;
            continue
        elseif (max(alt_array{p}) - min(alt_array{p})) < min_height
            prof_reject = bitset(prof_reject,6);
            db.reject{p} = prof_reject;
            continue
        end
        % Find all the pixels that impinge on this profile
        lat_logical = max(lattest) > min(corner_lat,[],1) & min(lattest) < max(corner_lat,[],1);
        lon_logical = max(lontest) > min(corner_lon,[],1) & min(lontest) < max(corner_lon,[],1);
        latlon_logical = lat_logical & lon_logical;
        omi_lon_p = omi_lon(latlon_logical); omi_lat_p = omi_lat(latlon_logical);
        loncorn_p = corner_lon(:,latlon_logical); latcorn_p = corner_lat(:, latlon_logical);
        omi_no2_p = omi_no2(latlon_logical); behr_no2_p = behr_no2(latlon_logical);
        tropopause_p = TropopausePres(latlon_logical); vza_p = vza(latlon_logical);
        terrain_pres_p = TerrainPres(latlon_logical);
        modis_cloud_p = modis_cloud(latlon_logical); total_omi_no2_p = total_omi_no2(latlon_logical);
        
        % Check each pixel for rejection criteria.
        pix_xx = true(size(omi_no2_p)); pix_reject = prof_reject*uint8(ones(size(pix_xx)));
        if isempty(pix_xx); pix_reject = bitset(prof_reject,4,1); end
        
        % Skip pixels with corners if the corners have differing signs and
        % are not near 0; this will confuse the algorithm that matches
        % spirals to pixels.  Check along the first dimension of the corner
        % variables.
        pix_xx = pix_xx & ~((any(abs(loncorn_p)>20) & abs(mean(sign(loncorn_p),1))~=1)');
        
        % Check the vza
        vza_xx = vza_p <= 60;
        pix_xx = pix_xx & vza_xx;
        pix_reject(~vza_xx) = bitset(pix_reject(~vza_xx),5*uint8(ones(size(find(~vza_xx)))),1);
        
        % Finally actually check if the profile falls in the pixel
        % using inpolygon(). Recall that we require there to be 20
        % valid measurements in the lowest 3 km.
        no2_3km = no2_array{p}(alt_array{p}<3);
        lon_3km = lon_array{p}(alt_array{p}<3);
        lat_3km = lat_array{p}(alt_array{p}<3);
        pix_coverage = zeros(size(omi_no2_p));
        for pix=1:size(loncorn_p,2)
            IN_3km = inpolygon(lon_3km, lat_3km, loncorn_p(:,pix), latcorn_p(:,pix));
            if sum(~isnan(no2_3km(IN_3km)))<numBLpoints
                pix_xx(pix) = false;
                pix_reject(pix) = bitset(pix_reject(pix),3,1);
                continue % Pixel must have 20 valid measurments between 0-3 km altitude (good sampling of boundary layer)
            end
            % Calculate what percentage of the profile actually falls in
            % this pixel; append to all values for this pixel
            IN_all = inpolygon(lon_array{p}, lat_array{p}, loncorn_p(:,pix), latcorn_p(:,pix));
            pix_coverage(pix) = sum(IN_all)/numel(no2_array{p});
        end
        
        % If no valid pixels are left, skip this profile.
        if sum(pix_xx)==0;
            continue
        end
        
        %Otherwise, cut down the vectors representing pixels that do match
        %this profile to the pixels that should be considered
        omi_lon_p = omi_lon_p(pix_xx); omi_lat_p = omi_lat_p(pix_xx);
        loncorn_p = loncorn_p(:,pix_xx); latcorn_p = latcorn_p(:, pix_xx);
        omi_no2_p = omi_no2_p(pix_xx); behr_no2_p = behr_no2_p(pix_xx);
        tropopause_p = tropopause_p(pix_xx); vza_p = vza_p(pix_xx);
        modis_cloud_p = modis_cloud_p(pix_xx); total_omi_no2_p = total_omi_no2_p(pix_xx);
        pix_coverage = pix_coverage(pix_xx); terrain_pres_p = terrain_pres_p(pix_xx);
        
        % Calculate the distance vector between the mean lat/lon of the
        % lowest 3 km of the profile and each pixel left.  This can be used
        % later for weighting.
        dist_vectors = zeros(size(omi_lat_p));
        prof_latlon = [nanmean(lon_3km), nanmean(lat_3km)];
        for a=1:numel(dist_vectors)
            omi_latlon = [omi_lon_p(a), omi_lat_p(a)]; 
            dist_vectors(a) = norm((prof_latlon - omi_latlon));
        end
        
        % Bin the NO2 data by pressure.
        [no2bins, presbins] = bin_omisp_pressure(pres_array{p}, no2_array{p});
        % Get the standard error of the bins
        [~,~,no2stderr] = bin_omisp_pressure(pres_array{p}, no2_array{p}, 'binmode','mean');
        % Bin the temperature
        [tempbins, temp_presbins] = bin_omisp_pressure(pres_array{p}, temp_array{p});
        
        
        % Get the top 10 NO2 measurements; if their median value is < 100
        % pptv, assume that the profile has sampled the free troposphere
        % and can be extrapolated to the tropopause safely.  If not, then
        % we will need to use the composite profile to fill in the bins
        % above the profile top. At the same time, get the bottom 10
        % NO2 measurements; we'll need them to extrapolate to the
        % surface.
        M = sortrows([pres_array{p}', no2_array{p}', radar_array{p}', temp_array{p}', alt_array{p}']);
        xx = find(~isnan(M(:,2)),10,'last'); zz = find(~isnan(M(:,2)),10,'first');
        bottom_med_no2 = median(M(xx,2)); top_med_no2 = median(M(zz,2));
        % Calculate the standard error for the top and bottom median points
        bottom_med_no2_stderr = std(M(xx,2))/sqrt(10);
        top_med_no2_stderr = std(M(zz,2))/sqrt(10);
        
        if DEBUG_LEVEL > 3;
            f2 = figure;
            plot(M(xx,2),M(xx,3),'color','b','linestyle','none','marker','o');
            line(M(zz,2),M(zz,3),'color','r','linestyle','none','marker','s');
            xlabel('NO2 mixing ratio'); ylabel('Radar altitude');
            title(sprintf('%s: Profile %d',Merge.metadata.date,profnum_array{p}));
            fprintf('\nPaused\n');
            pause;
            close(f2);
        end
        
        %Hains substitutes 1.5 ppt for any median no2 mixing ratios
        %less than the LoD (3 ppt)
        if top_med_no2 < 3;
            top_med_no2 = 1.5;
            if DEBUG_LEVEL>0; fprintf('Profile top median NO2 < LoD\n'); end
        end
        
        bottom_med_temp = nanmedian(M(xx,4)); top_med_temp = nanmedian(M(zz,4));
        bottom_med_radar_alt = nanmedian(M(xx,3)); bottom_med_pres = nanmedian(M(xx,1));
        bottom_med_GPS_alt = nanmedian(M(xx,5));
        % There is a chance that the radar system wasn't working at the
        % same time as the NO2 measurments, so if there were no
        % corresponding radar measurements in the lowest part of the
        % column, take whatever lowest 10 are available.
        if isnan(bottom_med_radar_alt);
            yy = find(~isnan(M(:,3)),10,'last');
            bottom_med_radar_alt = nanmedian(M(yy,3)); bottom_med_GPS_alt = nanmedian(M(yy,5));
            q_flag = bitset(q_flag,6,1);
        end
        % If the radar altitude is now a valid number, use it to get
        % the surface pressure.  Otherwise, use the GLOBE database and
        % find the nearest surface altitude, which must be converted to
        % kilometers
        if ~isnan(bottom_med_radar_alt)
            surface_alt = bottom_med_GPS_alt-bottom_med_radar_alt; surface_pres = 1013*exp(-surface_alt/7.4);
        else
            if DEBUG_LEVEL > 0; fprintf('  Retrieving GLOBE surface altitude.  May take a second...\n'); end
            surface_alt = nearest_GLOBE_alt(nanmean(lon_array{p}), nanmean(lat_array{p}))/1000;
            surface_pres = 1013*exp(-surface_alt/7.4);
            q_flag = bitset(q_flag,7,1);
        end
        
        if sum(~isnan(no2_array{p}))/numel(no2_array{p}) < 0.1
            q_flag = bitset(q_flag,4,1); % Set the 4th bit as a flag if less than 10% of the data points are non-fill values
        end
        
        if DEBUG_LEVEL > 3; f1 = figure; plot(no2bins, presbins, 'color','b','linewidth',10); end
        % If the top median no2 value is < 100 pptv, then it's safe to
        % assume that the profile has crossed the boundary layer and is
        % sampling the free troposphere.  If not, then extrapolating that
        % value to the tropopause will grossly overestimate the total
        % column.  In the latter case, append the composite profile on top
        % of the current one.
        if top_med_no2 < 100;
            no2bins(end) = top_med_no2;
            no2stderr(end) = top_med_no2_stderr;
            tempbins(end) = top_med_temp;
            topcol = [0 0.7 0];
        else
            % Get the altitude of the highest bin in the profile and find
            % all bins in the composite profile above that
            xx = pres_composite < min(presbins(~isnan(no2bins)));
            no2bins(xx) = no2_composite(xx);
            no2stderr(xx) = no2_composite_stderr(xx);
            tempbins(xx) = temp_composite(xx);
            topcol = [0.7 0 0.7];
            
            % Set the 3rd bit of the quality flag to 1 to indicate that a
            % composite column was appended
            q_flag = bitset(q_flag,3,1);
        end
        
        % Fill in any NaNs with interpolated values
        [tmp_pres, tmp_no2] = fill_nans(log(presbins),log(no2bins),'noclip');
        no2bins = exp(tmp_no2); presbins = exp(tmp_pres);
        [~, tempbins] = fill_nans(log(temp_presbins),tempbins,'noclip');
        
        if DEBUG_LEVEL > 3; figure(f1); line(no2bins,presbins,'color',topcol,'linewidth',4); end
        
        % Find the surface pressure, then compare it to the bottom bin with
        % NO2 measurements.  If it is above the bottom bin center, reset
        % that bin center to the surface pressure.  If below, then we'll
        % extrapolate the median lowest 10 NO2 measurements to surface
        % pressure.
        bb = find(~isnan(no2bins),1,'first');
        % Restrict the three bins to start from those that have NO2
        % values
        no2bins = no2bins(bb:end);
        presbins = presbins(bb:end);
        tempbins = tempbins(bb:end);
        no2stderr = no2stderr(bb:end);
        
        % If the surface pressure is less (i.e. above) the second
        % remaining bin, check with the user to proceed.
        if surface_pres < presbins(2);
            queststring = sprintf('Surface P (%.4f) less than second bin (%.4f). \nLow alt radar nan flag is %d. \n Continue?',surface_pres, presbins(2), bitget(q_flag,6));
            choice = questdlg(queststring,'Surface pressure','Yes','No','Abort run','No');
            switch choice
                case 'Yes'
                    cc = presbins < surface_pres;
                    no2bins = no2bins(cc);
                    presbins = presbins(cc);
                    tempbins = tempbins(cc);
                    no2stderr = no2stderr(cc);
                case 'No'
                    continue
                case 'Abort run'
                    error('spiral_ver:surface_pres','Surface pressure < second bin.');
            end
        elseif surface_pres < presbins(1); 
            % if surface is above the bottom bin center but below the
            % second bin, just reset the first bin pressure to be at the
            % surface.
            presbins(1) = surface_pres;
            nans = isnan(tempbins);
            tempbins(nans) = interp1(log(pres_composite(~nans)), tempbins(~nans), log(pres_composite(nans)),'linear','extrap'); % Just in case the temperature data doesn't cover all the remaining bins, interpolate it
        else
            % Otherwise, add a "surface bin" to interpolate to.
            no2nans = isnan(no2bins);
            no2bins = no2bins(~no2nans); tempbins = tempbins(~no2nans); presbins = presbins(~no2nans); no2stderr = no2stderr(~no2nans);
            no2bins = [bottom_med_no2, no2bins];
            no2stderr = [bottom_med_no2_stderr, no2stderr];
            tempbins = [bottom_med_temp, tempbins];
            presbins = [surface_pres, presbins];
        end
        
        
        if DEBUG_LEVEL > 3
            figure(f1); line(no2bins,presbins,'color','red');
            scatter_errorbars(no2bins,presbins,no2stderr,'direction','x','color','k','linewidth',2);
            title(sprintf('%s: profile %d',Merge.metadata.date,profnum_array{p}));
            set(gca,'YDir','reverse');
            fprintf('\nPaused\n');
            pause;
            close(f1);
        end
            
        % Insert the average OMI tropopause pressure of the pixels
        % considered as the final pressure bin center, then convert all
        % pressure bins to altitude, interpolate, and integrate.
        presbins(end) = nanmean(tropopause_p);
        altbins = -log(presbins ./ 1013) * 7.4;
        
        % Interpolate the NO2, temperature, and pressure data
        dz = 1; % integration segments in meters
        alt_profile = altbins(1):(dz/1000):altbins(end);
        no2_profile = interp1(altbins,no2bins,alt_profile,'linear');
        temp_profile = interp1(altbins,tempbins,alt_profile,'linear');
        pres_profile = exp(interp1(altbins,log(presbins),alt_profile,'linear')); % Linearly interpolate ln(P) since that is what depends linearly on altitude
        
        % Fill in the standard error nans, but we don't need to interpolate
        % to every 100 cm altitude point.
        [~,~,no2stderr] = fill_nans(altbins,no2bins,no2stderr,'noclip');
        
        % Carry out the numerical integration, 
        no2_column = 0;
        for z=1:numel(alt_profile)
            P_z = pres_profile(z); T = temp_profile(z); no2_z = no2_profile(z);
            conc_NO2 = (Av * P_z * no2_z * conv_fact)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
            no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
        end
        
        % The error propagation will be handled by considering this to be
        % effectively a trapezoid rule implementation.  The math behind
        % this is in my BEHR notebook from 4 Feb 2015 - Josh Laughner
        
        % Number density of air at each bin center
        Nair = (presbins .* Av)./(R .* tempbins);
        % Calculates the error for each trapezoid using vectors, then sum
        % up that vector and square root it to find the total column error.
        % Since the concentrations are defined in molec./cm^3, we need to
        % convert altitude bins from km to cm (hence the factor of 1e5).
        % Also, NO2 values are usually reported (by us at least) in pptv,
        % hence the conv_fact defaults to 1e-12.
        column_error = sqrt( sum( ((altbins(2:end) - altbins(1:end-1))*1e5/2).^2 .* (no2stderr(2:end).*conv_fact).^2 .* Nair(2:end).^2 +...
            ((altbins(2:end) - altbins(1:end-1))*1e5/2).^2 .* (no2stderr(1:end-1).*conv_fact).^2 .* Nair(1:end-1).^2 ) );
        % Double check that this number is a scalar. This will catch if a
        % matrix (rather than a vector) slips into this calculation.
        if ~isscalar(column_error)
            error(E.callError('''column_error'' is not a scalar value but it should be'));
        end
        
        
        % Save the output for this profile
        prof_lon_out(p) = nanmean(lon_array{p});
        prof_lat_out(p) = nanmean(lat_array{p});
        omi_no2_out(p) = nanmean(omi_no2_p);
        behr_no2_out(p) = nanmean(behr_no2_p);
        air_no2_out(p) = no2_column;
        
        % Bin the aerosol data for this profile and find the max extinction
        % value and total integrated aerosol extinction.
        if aerfield ~= 0;
            % Check the quality of the aerosol data in this profile. If the
            % entire profile is NaNs, skip it
            aer_quality = uint8(0);
            if all(isnan(aer_array{p}))
                % Set the first bit in aer_quality to 1 if all data is NaNs
                aer_quality = bitset(aer_quality,1);
                aer_max = -9e9;
                aer_int = -9e9;
            else
                if sum(isnan(aer_array{p}))/numel(aer_array{p}) > 0.9
                    % Set the second bit if >90% of the profile has NaNs
                    % for aerosol data.
                    aer_quality = bitset(aer_quality,2);
                end
                [aerbins] = bin_omisp_pressure(pres_array{p}, aer_array{p});
                aer_max = max(aerbins);
                if isempty(aer_max);
                    aer_max = -9e9;
                end
                
                binwidth = 0.25; % km
                [aerintbins, aerbinmid] = bin_vertical_profile(alt_array{p}, aer_array{p}, binwidth);
                [aerbinmid, aerintbins] = fill_nans(aerbinmid,aerintbins);
                if numel(aerintbins) == 1
                    % Set the third bit if only one bin had actual data
                    aer_quality = bitset(aer_quality,3);
                end
                % Integrate, converting the aerosol extinction from Mm^-1
                % to m^-1 and the bin altitudes from km to m.
                aer_int = trapz(aerbinmid*1e3, aerintbins*1e-6);
            end
        end
        
        % If any bits in the quality flag are set, set the summary
        % bit; then append the quality flag to all those for this pixel
        if any(q_flag); q_flag = bitset(q_flag,1,1); end
        
        db.all_omi{p} = omi_no2_p;
        db.all_behr{p} = behr_no2_p;
        db.quality_flags{p} = q_flag;
        db.coverage_fraction{p} = pix_coverage;
        db.dist_vectors{p} = dist_vectors;
        db.latcorn{p} = latcorn_p;
        db.loncorn{p} = loncorn_p;
        db.strat_NO2{p} = total_omi_no2_p - omi_no2_p;
        db.modis_cloud{p} = modis_cloud_p;
        db.profnums{p} = profnum_array{p};
        db.reject{p} = pix_reject;
        db.lon_3km{p} = lon_3km;
        db.lat_3km{p} = lat_3km;
        if aerfield ~= 0; 
            db.aer_max_out{p} = aer_max; 
            db.aer_int_out{p} = aer_int;
            db.aer_quality{p} = aer_quality;
        end
        db.column_error{p} = column_error;
        
    end % End the loop over all profiles

    % Clean up the output variables
    if clean_bool
        fills = prof_lon_out == -9e9;
        prof_lon_out = prof_lon_out(~fills);
        prof_lat_out = prof_lat_out(~fills);
        omi_no2_out = omi_no2_out(~fills);
        behr_no2_out = behr_no2_out(~fills);
        air_no2_out = air_no2_out(~fills);
        
        db.all_omi = db.all_omi(~fills);
        db.all_behr = db.all_behr(~fills);
        db.quality_flags = db.quality_flags(~fills);
        db.coverage_fraction = db.coverage_fraction(~fills);
        db.dist_vectors = db.dist_vectors(~fills);
        db.latcorn = db.latcorn(~fills);
        db.loncorn = db.loncorn(~fills);
        db.strat_NO2 = db.strat_NO2(~fills);
        db.modis_cloud = db.modis_cloud(~fills);
        db.profnums = db.profnums(~fills);
        db.reject = db.reject(~fills);
        db.lon_3km = db.lon_3km(~fills);
        db.lat_3km = db.lat_3km(~fills);
        if aerfield ~= 0; 
            db.aer_max_out = db.aer_max_out(~fills); 
            db.aer_int_out = db.aer_int_out(~fills);
            db.aer_quality = db.aer_quality(~fills);
        end
        db.column_error = db.column_error(~fills);
    end
end
end
