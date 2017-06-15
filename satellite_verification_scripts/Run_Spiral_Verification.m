% Run_Spiral_Verification
%
%   This script will automatically execute sprial_verification.m for each
%   swath of each day within the range of dates given.
%
%   RUN_SPIRAL_VERIFICATION() - uses the default campaign_name and will
%   include all profiles. Output is placed directly in base workspace.
%
%   RUN_SPIRAL_VERIFICATION('all', CAMPAIGN_NAME) - includes all profiles
%   for the given campaign. See MERGE_FIELD_NAMES for valid campaign names.
%   This format is necessary to maintain compatibility with
%   RUN_ALL_AER_CATEGORIES.
%
%   RUN_SPRIAL_VERIFICATION(PROFNUMS) - only includes the specified
%   profiles for the default campaign. Output is placed directly in the
%   base workspace.
%
%   RUN_SPIRAL_VERIFICATION(PROFNUMS, CAMPAIGN_NAME) - only includes the
%   specified profiles for the given campaign. See MERGE_FIELD_NAMES for
%   valid campaign names. Output is placed directly in base workspace.
%
%   [LON, LAT, OMINO2, BEHRNO2, AIRNO2, DB, DATES] =
%   RUN_SPIRAL_VERIFICATION(...) with any previous syntax, outputs the
%   variables as a normal function. LON, LAT are the average profile
%   lon/lat, OMINO2 and BEHRNO2 are the matched satellite NO2 columns for
%   NASA SP and BEHR respectively. AIRNO2 is the aircraft VCDs calculated
%   from the profiles. DB is a structure containing much extra information,
%   and DATES is a cell array of the dates for each profiles.
%
%   In any case, only profiles that intersect at least one valid OMI pixel
%   are used. 
%
%  Josh Laughner <joshlaugh5@gmail.com> 4 Jul 2014

function varargout = Run_Spiral_Verification(profnums_in, campaign_name)

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% The campaign name. Current valid strings are 'discover-md','discover-ca',
% 'discover-tx', 'seac4rs', 'dc3', 'arctas-b', 'arctas-carb'. This is used
% to automatically find the campaign dates, the campaign directory, and the
% data field names. If you don't want to retrieve this automatically, set
% this to an empty string.
if ~exist('campaign_name','var')
    campaign_name = 'discover-md'; % Which campaign this is for. Used to automatically find field names
end

% Grab the dates and directory for the campaign unless the campaign name is
% empty.
if ~isempty(campaign_name)
    [~, dates, directory, range_file] = merge_field_names(campaign_name);
end

% The dates to run - leave as empty strings to automatically do the entire
% campaign, which is read from the merge_field_names function output above.
date_start = '';
date_end = '';

if isempty(date_start) || isempty(date_end)
    fprintf('Setting start and end dates based on the campaign name.\n');
    date_start = dates{1};
    date_end = dates{2};
end

% The directories where to find the data needed for the run.  The merge
% directory (where the campaign data is located) can be automatically
% determined from the campaign name (leave the string empty if this is
% desired), but for now the BEHR directory will be set manually.  The BEHR
% prefix is the part of the BEHR file name before the date; wildcards are
% allowed so long as there is enough of the prefix there to uniquely
% identify 1 file per date in the given directory.
merge_dir = '';
behr_dir = '/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF';
behr_prefix = 'OMI_BEHR_v2-1B_';


if isempty(merge_dir)
    fprintf('Setting merge file directory based on the campaign name.\n');
    merge_dir = directory;
end

% Start and end times (in military format) for which profiles to consider.
% General recommendation is +/-1.5 hr from overpass.
starttime = '12:00';
endtime = '15:00';

% Time zone (3 letter abbreviation). Set to 'auto' to determine based on
% the longitude of the data

tz = 'auto';


% Which fields from the merge files to use. Set them to empty strings ('')
% to automatically guess the correct field for the given campaign.  

no2field = 'lif'; % Which NO2 data field to use. Leave as empty string or 'lif' for our LIF data, set to 'cl' for chemiluminescence data, or any other string to override.
conv_fact = 1e-12; % Conversion factor for NO2 data from part-per-whatever to part-per-part. Usually 1e-12, i.e. NO2 data is in pptv.
aerfield = ''; % Which aerosol extinction field to use.
ssafield = ''; % Which aerosol SSA field to use.
altfield = ''; % Which altitude field to use. Can set to 'pressure' or 'gps' ('' defaults to gps), or override.
radarfield = ''; % The field for radar altitude.

% Variables to allow or disallow the use of a profile
min_height = 0; % The minimum difference between the top and bottom of a profile. Set to 0 to ignore (max - min > 0 always).
numBLpoints = 20; % The number of data points required in the bottom 3 km to ensure good BL sampling. Hains et. al. recommends 20.
minRadarAlt = 0.5; % Height above the surface (in km) a profile must be below to ensure good BL sampling. Hains et. al. recommends 0.5 km (500 m).

% Set to 1 to include ground site data, or 0 to use only aircraft data.
useground = 1;

% These variables are used to subset the aircraft data into profiles used
% to generate column data to compare against satellite data. 
%   'profiles' should either be set to: 
%       1) an empty string or the name of the profile number field 
%       2) the string 'ranges'. 
%   Option 1 only works with the DISCOVER campaigns so far, as those are the
% only campaigns that have specifically generated spirals labeled with
% profile numbers.  An empty string will automatically detect the correct
% field name for a campaign.  Using this option, the "profnums" variable
% can be set with an empty matrix to use all profiles, or to a matrix
% containing specific profile numbers to be examined. Such a matrix can be
% generated using categoriz_aerosol_profile in the Aerosol Effects/Profile
% Classification directory.
%   Option 2 requires a range file (a file containing the "Ranges"
% structure generated by "record_profile_ranges" in the utility_scripts
% folder) and for the path to that file to be specified.  This file
% contains a list of UTC ranges that have been specified as profiles of
% interest.  If range_file is left blank, the program will look at the
% range files returned from merge_field_names.  If there is one, that one
% will be used, otherwise the user is presented with his options.
profile_input = '';
profnums = 'fetch'; % set to 'fetch' to use the first input to this function as the profile numbers


% Fields for pressure and temperature, set to empty strings to
% automatically detect.
presfield = '';
tempfield = '';

% Satellite variables
cloud_product = 'omi'; % Can be 'omi', 'modis', or 'rad'
cloud_frac_max = 0.2; % Maximum cloud fraction to allow in a pixel. Recommended 0.2 for OMI, 0 for MODIS, and 0.5 for radiance.
row_anomaly = 'XTrackFlags'; % How to reject for the row anomaly - can be 'AlwaysByRow', 'RowsByTime', 'XTrackFlags', and 'XTrackFlagsLight'
behrfield = 'BEHRColumnAmountNO2Trop'; % The field in the Data structure with NO2 column data. 
                                          % Choices include 'ColumnAmountNO2Trop' (OMI SP column), 'BEHRColumnAmountNO2Trop' (BEHR column) 
                                          % and 'BEHR_R_ColumnAmountNO2Trop' only available in files where the column was reprocessed with
                                          % an AMF derived from in-situ measurements.

useghost = 0;   % set to 0 to use visible columns only (in new BEHR)
                % set to 1 to multiply columns by the ghost factor (if available)
                % set to 2 to divide columns by g.f. (old method)
% Debugging variables
DEBUG_LEVEL = 2; % This will also be passed to the spiral verification function.
clean = 1; % Set to 0 to keep all pixel comparisons, even those with fill values. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA PREPARATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If we are using ranges, load the range file and extract the date
% information - we'll need that to match up the correct set of ranges with
% the correct date.    
fprintf('Using %s as the range file\n',range_file);


if strcmpi(profile_input, 'ranges')
    load(range_file); % loads the Ranges variable, which is a data structure
    range_dates = cellstr(datestr({Ranges.Date},29));
    range_bool = true;
else
    range_bool = false;
end

% If the user wants to use profnums input to the function, retrieve them
% here. If there is no input, pass an empty matrix.
if strcmpi(profnums,'fetch') && nargin > 0
    if isempty(profnums_in)
        % Use a random number here that the profile number will never be
        % because if an empty set of profile numbers is given, it indicates
        % that no profiles should be matched; but spiral_verf. uses an
        % empty matrix to indicate that it shouldn't filter the profile
        % numbers at all.
        if strcmpi(profile_input, 'ranges')
            profnums = [-127, -127];
        else
            profnums = -127;
        end
    elseif strcmpi(profnums_in,'all')
        profnums = [];
    else
        profnums = profnums_in;
    end
elseif nargin > 0 && ~strcmp(profnums,'fetch')
    warning('Input detected, but profnums is not set to fetch that input.');
elseif strcmpi(profnums, 'fetch')
    profnums = [];
end

% This boolean will help identify whether the run succeeded, that is, if
% any files where actually loaded.
run_bool = false;


%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

dates = datenum(date_start):datenum(date_end);

S=0; %clear('dbs');
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('%s%s%s%s.mat',behr_prefix,year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    behr_files = dir(fullfile(behr_dir,behr_filename));
    if numel(behr_files)==1
        load(fullfile(behr_dir,behr_files(1).name),'Data')
    elseif isempty(behr_files)
        if DEBUG_LEVEL > 1; fprintf('No BEHR file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of BEHR files for %s is not 1 or 0',datestr(dates(d)));
    end

    % Find the UTC range data for this date, if we're using ranges instead
    % of profile numbers
    if range_bool
        xx = find(strcmp(curr_date,range_dates));
        if isempty(xx);
            if DEBUG_LEVEL > 0; fprintf('No UTC ranges found for %s, skipping\n',curr_date); end
            continue
        end
        profile_input = Ranges(xx).Ranges;
        
        % As we start considering very specific events (e.g. fires) not all
        % days will have ranges defined. 
        if isempty(profile_input); continue; end
    end
    
    % If no BEHR or Merge files are ever loaded, we won't get to this line,
    % so this will be false. This boolean will be used to call an error an
    % little later if still false, and that error will avoid the more
    % confusing error that dbs doesn't exist.
    if ~run_bool; run_bool = true; end
    
    for swath=1:numel(Data)
        S=S+1;
        [lon_i{S}, lat_i{S}, omino2_i{S}, behrno2_i{S}, airno2_i{S}, dbs(S)] = spiral_verification_avg_pix2prof(Merge,Data(swath),tz,...
            'behrfield',behrfield,...
            'starttime',starttime,... 
            'endtime',endtime,... 
            'profiles',profile_input,...
            'profnums',profnums,... 
            'campaign_name',campaign_name,...
            'no2field',no2field,... 
            'conv_fact',conv_fact,...
            'aerfield',aerfield,... 
            'ssafield',ssafield,...
            'altfield',altfield,... 
            'radarfield',radarfield,... 
            'presfield',presfield,... 
            'tempfield',tempfield,... 
            'cloud_product',cloud_product,... 
            'cloud_frac_max',cloud_frac_max,... 
            'rowanomaly',row_anomaly,... 
            'min_height',min_height,...
            'numBLpoints',numBLpoints,...
            'minRadarAlt',minRadarAlt,...
            'useground',useground,...
            'useghost',useghost,...
            'DEBUG_LEVEL',DEBUG_LEVEL,... 
            'clean',clean); 
        date_cell{S} = repmat({curr_date},numel(lon_i{S}),1);
    end
end

if ~run_bool;
    error(E.callError('run_failure','The loop never executed completely, meaning in no case was both a BEHR and Merge file loaded. Check the dates, file paths, and BEHR prefix.'));
end

% concatenate the output
[db_iall, lon_iall, lat_iall, omino2_iall, behrno2_iall, airno2_iall, dates_iall] = match_arrays2db(dbs, lon_i, lat_i, omino2_i, behrno2_i, airno2_i, date_cell);

% Save the inputs to the spiral_verification function; this way by saving
% the db_iall structure we can recreate this run later if necessary.
db_iall.run.date_start = date_start;
db_iall.run.date_end = date_end;
db_iall.run.merge_dir = merge_dir;
db_iall.run.behr_dir = behr_dir;
db_iall.run.behr_prefix = behr_prefix;
db_iall.run.behrfield = behrfield;
db_iall.run.range_file = range_file;
db_iall.run.starttime = starttime;
db_iall.run.endtime = endtime;
db_iall.run.timezone = tz;
if range_bool;
    db_iall.run.profile_input = 'ranges';
else
    db_iall.run.profile_input = profile_input;
end
db_iall.run.profnums = profnums;
db_iall.run.campaign_name = campaign_name;
db_iall.run.no2field = no2field;
db_iall.run.aerfield = aerfield;
db_iall.run.altfield = altfield;
db_iall.run.radarfield = radarfield;
db_iall.run.presfield = presfield;
db_iall.run.tempfield = tempfield;
db_iall.run.cloud_product = cloud_product;
db_iall.run.cloud_frac_max = cloud_frac_max;
db_iall.run.row_anomaly = row_anomaly;
db_iall.run.min_height = min_height;
db_iall.run.numBLpoints = numBLpoints;
db_iall.run.minRadarAlt = minRadarAlt;
db_iall.run.useground = useground;
db_iall.run.useghost = useghost;


if nargout == 0;
    putvar(db_iall, lon_iall, lat_iall, omino2_iall, behrno2_iall, airno2_iall, dates_iall);
else
    varargout{1} = lon_iall;
    varargout{2} = lat_iall;
    varargout{3} = omino2_iall;
    varargout{4} = behrno2_iall;
    varargout{5} = airno2_iall;
    varargout{6} = db_iall;
    varargout{7} = dates_iall;
end
end
