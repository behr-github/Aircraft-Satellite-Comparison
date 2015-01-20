% Run_Spiral_Verification
%
%   This script will automatically execute sprial_verification.m for each
%   swath of each day within the range of dates given.
%
%  Josh Laughner <joshlaugh5@gmail.com> 4 Jul 2014

function varargout = Run_Spiral_Verification(profnums)

date_start = '07/01/2011';
date_end = '07/31/2011';

no2field = 'NO2_LIF';
altfield = 'GPS_ALT';

starttime = '12:00';
endtime = '15:00';

tz = 'auto';
%profnums = profstr.NO2AboveLow;
%profnums = [2012, 6011, 3026, 5027, 5028, 6024, 6027, 3042, 6038]; %coincident
%profnums = [2021, 2022,3023,4021,6021,8003,2023,2036,3040,4034,5041]; %shielding

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
%behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR_REPROCESSED/';
behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR/';
behr_prefix = 'OMI_BEHR_omi*';

DEBUG_LEVEL = 0;

if nargin < 1; profnums = 1:1e7; end

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

    
    for swath=1:numel(Data)
        S=S+1;
        [lon_i{S}, lat_i{S}, omino2_i{S}, behrno2_i{S}, airno2_i{S}, dbs(S)] = spiral_verification_avg_pix2prof(Merge,Data(swath),tz,'DEBUG_LEVEL',DEBUG_LEVEL,...
            'no2field',no2field,'altfield',altfield,'starttime',starttime,'endtime',endtime,'rowanomaly','XTrackFlags','behrfield','BEHRColumnAmountNO2Trop','profnums',profnums);
        date_cell{S} = repmat({curr_date},numel(lon_i{S}),1);
    end
end

% concatenate the output
[db_iall, lon_iall, lat_iall, omino2_iall, behrno2_iall, airno2_iall, dates_iall] = match_arrays2db(dbs, lon_i, lat_i, omino2_i, behrno2_i, airno2_i, date_cell);
% lon_iall = cat(1,lon_i{:});
% lat_iall = cat(1,lat_i{:});
% omino2_iall = cat(1, omino2_i{:});
% behrno2_iall = cat(1, behrno2_i{:});
% airno2_iall = cat(1, airno2_i{:});
% dates_iall = cat(1,date_cell{:});
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