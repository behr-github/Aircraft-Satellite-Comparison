% Run_Spiral_Verification
%
%   This script will automatically execute sprial_verification.m for each
%   swath of each day within the range of dates given.  This version of the
%   script will pass the UTC 
%
%  Josh Laughner <joshlaugh5@gmail.com> 4 Jul 2014

date_start = '03/01/2006';
date_end = '05/15/2006';

debugreject = 0; %Set to 1 to enable looking at pixel rejection.

no2field = 'NO2';
altfield = 'ALTITUDE_GPS';
radarfield = 'ALTITUDE_RADAR';
tempfield = 'TEMP_STAT_C';
presfield = 'STAT_PRESSURE';

tz = 'auto'; %set to 'auto' to calculate the time zone based on the mean longitude of the flight

merge_dir = '/Volumes/share/GROUP/INTEX-B/Matlab files/';
behr_dir = '/Volumes/share-sat/SAT/OMI/Bare_SP_Files/';
range_file = '/Volumes/share/GROUP/INTEX-B/INTEXB_Profile_UTC_Ranges_Inclusive.mat';

DEBUG_LEVEL = 1;

load(range_file); range_dates = cellstr(datestr({Ranges.Date},29));
dates = datenum(date_start):datenum(date_end);

S=0; clear('db'); date_list = cell(numel(dates)*15);
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('OMI_SP_*%s%s%s.mat',year,month,day);
    
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
    
    % Find the UTC range data for this date
    xx = find(strcmp(curr_date,range_dates));
    if isempty(xx);
        error('run_spiral:ranges','No UTC ranges found for %s',curr_date);
    end
    
    for swath=1:numel(Data)
        S=S+1;
        [lon_radrowanyvza{S}, lat_radrowanyvza{S}, omino2_radrowanyvza{S}, behrno2_radrowanyvza{S}, airno2_radrowanyvza{S}, db(S)] = spiral_verification_avg_pix2prof(Merge,Data(swath),tz,'DEBUG_LEVEL',1,'no2field',no2field,'profiles',Ranges(xx).Ranges,'cloud_product','rad','cloud_frac',0.5,'radarfield',radarfield,'altfield',altfield,'presfield',presfield,'tempfield',tempfield,'rowanomaly','RowsByTime','clean',~debugreject);
        date_list{S} = curr_date;
    end
    dummy = 1;
end

date_list = date_list(1:S);

% concatenate the output
lon_radrowanyvzaall = cat(1,lon_radrowanyvza{:});
lat_radrowanyvzaall = cat(1,lat_radrowanyvza{:});
omino2_radrowanyvzaall = cat(1, omino2_radrowanyvza{:});
behrno2_radrowanyvzaall = cat(1, behrno2_radrowanyvza{:});
airno2_radrowanyvzaall = cat(1, airno2_radrowanyvza{:});