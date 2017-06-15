function read_merge_data
% read_merge_data: reads in any merge files from flight campaigns and saves
% a .mat file containing a data structure for a given flight.  Said data
% structure will contain the each data set as a sub field, with fill value
% and units.
%
%   Josh Laughner <joshlaugh5@gmail.com> 22 May 2014

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  USER INPUT   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Set this to 1 to select icart and save directories using a standard file
% browser window

user_select_dir = 1;

% The location of the campaign (Baltimore_DC, Texas, CA, Colorado). Must
% match directory structure
location = 'Texas';

% The file types (e.g. 1 sec merge). Must match directory structure.
data_type = '1 sec merges';

% The name of the overall directory containing the files to be read in.
% Location and merge type will be populated
icart_path = '/Volumes/share/GROUP/DISCOVER-AQ/';

% The directory to save the matlab files to
save_path = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft';

% Set this to 1 to run a single file, set to 0 to run a full directory
single_file = 0;
% The date of the file to run, only has an effect if single file is set to 1
filedate = '07/01/2011';

% Level of output to console; 0 = nothing, 1 = minimal, 2 = all messages
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if user_select_dir;
    icart_dir = uigetdir('/Volumes','Select the directory with the ICART files in it');
    fprintf('Will read icart files from %s\n',icart_dir);
    save_path = uigetdir('/Volumes','Select the directory to save the resulting .mat files to');
    fprintf('Will save mat files to %s\n',save_path);
    save_bool = 1;
    if all(icart_dir == 0) || all(save_path == 0); error('read_merge_data:user_cancel','User canceled run.'); end
    
    savetitle = sprintf('Enter the save file name. This will preceed the data date. Canceling will use %s',[location,'_',data_type]);
    user_savename = inputdlg(savetitle); user_savename = user_savename{1};
    if ~isempty(user_savename) && strcmp(user_savename(end),'_');
        user_savename = user_savename(1:end-1);
    end
else % If using the coded directory, validate said directories
    % Check that the save path exists. Ask the user if they want to create the
    % folder, abort, or don't save.
    if ~exist(save_path,'dir')
        response = input('Save folder does not exist. [C]reate, [D]on''t save, or abort (default):  ', 's');
        if strcmpi(response,'c'); mkdir(save_path); save_bool = 1;
        elseif strcmpi(response,'d'); save_bool = 0;
        else error('read_merge:save','User aborted: save directory does not exist');
        end
    else
        save_bool = 1;
    end
    
    icart_dir = fullfile(icart_path, location, data_type);
    % Check that the data directory exists
    if ~exist(icart_dir,'dir');
        error('read_icart:data_dir_DNE','Data directory does not exist.');
    end
end



% Collect the list of files to read; if single_file is 1, restrict them to
% just that file.  In either case, restrict it to .ict files (this avoids
% getting hidden/unwanted files)

if single_file
    filedatestr = datestr(filedate,29);
    fileyear = filedatestr(1:4);
    filemonth = filedatestr(6:7);
    fileday = filedatestr(9:10);
    filename = ['*',fileyear,filemonth,fileday,'*.ict'];
    merge_files = dir(fullfile(icart_dir, filename));
else
    merge_files = dir(fullfile(icart_dir, '*.ict'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   MAIN LOOP   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through all files in the folder
numfiles = numel(merge_files);

for a = 1:numfiles
    if DEBUG_LEVEL > 0; fprintf('Reading in file %s\n', merge_files(a).name); end
    
    [Merge, header, DataTable] = read_icartt_file(fullfile(icart_dir, merge_files(a).name));
    
    % Save the day's structure, table, and header as a .mat file
    if save_bool
        file_date = regexprep(curr_date,'/','');
        if user_select_dir && ~isempty(user_savename);
            savename = [user_savename,'_',file_date];
        else
            savename = [location, '_', data_type, '_', file_date];
        end
        savename = regexprep(savename,'\W','_');
        
        savenum = 2;
        while true
           if ~exist(fullfile(save_path,savename),'file'); break; end
           savename = [savenum, '_',num2str(savenum)];
           savenum = savenum + 1;
        end
        
        save(fullfile(save_path,savename),'Merge','header','DataTable');
        if DEBUG_LEVEL > 0; fprintf('    File saved as %s in %s\n',savename,save_path); end
    end
end
end
