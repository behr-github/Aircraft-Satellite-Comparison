% read_merge_data: reads in any merge files from flight campaigns and saves
% a .mat file containing a data structure for a given flight.  Said data
% structure will contain the each data set as a sub field, with fill value
% and units.
%
%   Josh Laughner <joshlaugh5@gmail.com> 22 May 2014

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  USER INPUT   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% The location of the campaign (Baltimore_DC, Texas, CA, Colorado). Must
% match directory structure
location = 'Baltimore_DC';

% The file types (e.g. 1 sec merge). Must match directory structure.
data_type = 'Ozone sondes';

% The name of the overall directory containing the files to be read in.
% Location and merge type will be populated
icart_path = '/Volumes/share/GROUP/DISCOVER-AQ/';

% The directory to save the matlab files to
save_path = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Sondes';

% Set this to 1 to run a single file, set to 0 to run a full directory
single_file = 0;
% The date of the file to run, only has an effect if single file is set to 1
filedate = '07/01/2011';

% Level of output to console; 0 = nothing, 1 = minimal, 2 = all messages
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Check that the data directory exists
if ~exist(fullfile(icart_path, location, data_type),'dir');
    error('read_icart:data_dir_DNE','Data directory does not exist.');
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
    merge_files = dir(fullfile(icart_path, location, data_type, filename));
else
    merge_files = dir(fullfile(icart_path, location, data_type, '*.ict'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   MAIN LOOP   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through all files in the folder
numfiles = numel(merge_files);
for a = 1:numfiles
    if DEBUG_LEVEL > 0; fprintf('Reading in file %s\n', merge_files(a).name); end
    
    % Create the data structure
    Merge = struct;
    
    Merge.metadata.file = merge_files(a).name;
    
    fid = fopen(fullfile(icart_path, location, data_type, merge_files(a).name));
    N = fscanf(fid, '%d, %d\n'); % Get the number of lines before the data starts
    
    line_count = 2; % Keep track of what line we're on
    
    % Skip lines 2 & 3
    for l=2:3; fgetl(fid); line_count = line_count + 1; end
    
    % Get the description of the file for the metadata
    Merge.metadata.description = fgetl(fid); line_count = line_count + 1;
    % Create the UTC unit field
    Merge.Data.UTC.fill = 'N/A';
    Merge.Data.UTC.unit = 'seconds after midnight in UTC time';
    
    % Skip lines 5-6
    for l=5:6; fgetl(fid); line_count = line_count + 1; end
    
    % Get the data date
    y = fscanf(fid,'%d, %d, %d');
    curr_date = datestr([num2str(y(1)),'/',num2str(y(2)),'/',num2str(y(3))],29);
    Merge.metadata.date = curr_date;
    
    % Skip the rest of line 7-9
    for l=7:9; fgetl(fid); line_count = line_count + 1; end
    
    % Get the number of fields we'll need to handle
    numfields = str2double(fgetl(fid)); line_count = line_count + 1;
    
    % Skip one line
    fgetl(fid); line_count = line_count + 1;
    
    % Read in the fill values
    fill_vals = fscanf(fid,[repmat('%d, ',1,numfields-1),'%d\n']);
    line_count = line_count + 1;
    
    % For each field, read the unit and fill value into the data structure
    for f = 1:numfields
        line = fgetl(fid); line_count = line_count + 1;
        comma_pos = strfind(line,',');
        field = line(1:(comma_pos-1));
        unit = line((comma_pos+2):end);
        
        % Sanitize the field name, characters that don't belong in variable
        % names will be removed
        field = regexprep(field,'\W','');
        eval(sprintf('Merge.Data.%s.Unit = ''%s'';',field,unit))
        eval(sprintf('Merge.Data.%s.Fill = %d;',field,fill_vals(f)));
    end
    
    % Read lines until we've hit the big table of values
    while line_count < N(1)
        line = fgetl(fid);
        %Parse the line to see if it has any useful data
        if length(line) < 12;
        elseif strcmp(line(1:9),'DATA_INFO'); Merge.metadata.info = line(11:end);
        elseif strcmp(line(1:9),'ULOD_FLAG'); Merge.metadata.upper_lod_flag = str2double(line(11:end));
        elseif strcmp(line(1:10), 'ULOD_VALUE'); Merge.metadata.upper_lod_value = line(12:end);
        elseif strcmp(line(1:9), 'LLOD_FLAG'); Merge.metadata.lower_lod_flag = str2double(line(11:end));
        elseif strcmp(line(1:10), 'LLOD_VALUE'); Merge.metadata.lower_lod_value = line(12:end);
        end
        
        line_count = line_count + 1;
    end
    
    % Read in the header row
    header = cell(1,numfields+1);
    for f = 1:numfields+1
        field = fscanf(fid,'%s,');
        field = regexprep(field,'\W','');
        header{f} = field;
    end
    
    % Close the fid; we will use dlmread now and no longer need this file
    % open
    fclose(fid);
    
    if DEBUG_LEVEL > 0; fprintf('   Loading the data table\n'); end
    % Read in the full table
    DataTable = dlmread(fullfile(icart_path, location, data_type, merge_files(a).name),',',N(1),0);
    
    % Double check that the table size matches the header row
    if size(DataTable,2) ~= numel(header);
        error('read_merge:data_import','Table imported and header row do not have the same number of columns');
    end
    
    % Append the table data to the Merge structure
    for r = 1:numfields+1
        field = header{r};
        if DEBUG_LEVEL > 1; fprintf('\tField: %s', field); end
        data_column = DataTable(:,r)';
        eval(sprintf('Merge.Data.%s.Values = data_column;',field));
        if DEBUG_LEVEL > 1; fprintf(' Done.\n'); end
    end
    
    % Save the day's structure, table, and header as a .mat file
    if save_bool
        file_date = regexprep(curr_date,'/','');
        savename = [location, '_', data_type, '_', file_date];
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
