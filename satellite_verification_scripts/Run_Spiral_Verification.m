% Run_Spiral_Verification
%
%   This script will automatically execute sprial_verification.m for each
%   swath of each day within the range of dates given.
%
%  Josh Laughner <joshlaugh5@gmail.com> 4 Jul 2014

date_start = '07/01/2011';
date_end = '07/31/2011';

no2field = 'NO2_LIF';

starttime = '13:00';
endtime = '14:00';

tz = 'auto';

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR/';

DEBUG_LEVEL = 0;

dates = datenum(date_start):datenum(date_end);

S=0; clear('db');
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('OMI_BEHR_*%s%s%s.mat',year,month,day);
    
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

    
    for swath=1:numel(Data)
        S=S+1;
        [lon_i{S}, lat_i{S}, omino2_i{S}, behrno2_i{S}, airno2_i{S}, db(S)] = spiral_verification(Merge,Data(swath),tz,'DEBUG_LEVEL',DEBUG_LEVEL,'no2field',no2field,'starttime',starttime,'endtime',endtime);
        date_cell{S} = repmat({curr_date},numel(lon_i{S}),1);
    end
end

% concatenate the output
lon_iall = cat(1,lon_i{:});
lat_iall = cat(1,lat_i{:});
omino2_iall = cat(1, omino2_i{:});
behrno2_iall = cat(1, behrno2_i{:});
airno2_iall = cat(1, airno2_i{:});
dates_iall = cat(1,date_cell{:});