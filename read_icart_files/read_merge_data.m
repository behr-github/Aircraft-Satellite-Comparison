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
    
    % Create the data structure
    Merge = struct;
    
    Merge.metadata.file = merge_files(a).name;
    
    fid = fopen(fullfile(icart_dir, merge_files(a).name));
    line1 = fgetl(fid); % Some merge files separate the two numbers on the first line with a comma, some don't. This gets the whole line...
    N = sscanf(line1, '%d'); % ...then this gets the first number on that line.
    
    if length(N) == 1; delim = ',';
    else delim = '';
    end
    
    line_count = 2; % Keep track of what line we're on
    
    % Skip lines 2 & 3
    for l=2:3; fgetl(fid); line_count = line_count + 1; end
    
    % Get the description of the file for the metadata
    Merge.metadata.description = fgetl(fid); line_count = line_count + 1;
    % Create the UTC unit field
    Merge.Data.UTC.fill = 'N/A';
    Merge.Data.UTC.Unit = 'seconds after midnight in UTC time';
    
    % Skip lines 5-6
    for l=5:6; fgetl(fid); line_count = line_count + 1; end
    
    % Get the data date.  Check that y contains at least 3 values, if it
    % doesn't, the file probably uses space delimiters instead of commas.
    line7 = fgetl(fid); line_count = line_count + 1;
    y = sscanf(line7,'%d, %d, %d');
    if length(y) < 3
        y = sscanf(line7,'%d %d %d');
    end
    curr_date = datestr([num2str(y(1)),'/',num2str(y(2)),'/',num2str(y(3))],29);
    Merge.metadata.date = curr_date;
    
    % Skip lines 8-9
    for l=8:9; fgetl(fid); line_count = line_count + 1; end
    
    % Get the number of fields we'll need to handle
    numfields = str2double(fgetl(fid)); line_count = line_count + 1;
    if a==1; field_bool = -1*ones(1,numfields+1); end
    % Skip one line
    fgetl(fid); line_count = line_count + 1;
    
    % Read in the fill values, creating the format spec appropriately based
    % on the delimiter decided from the first line.
    line12 = fgetl(fid);
    if strcmpi(delim,',');
        fill_vals = sscanf(line12,repmat('%d, ',1,numfields));
    else
        fill_vals = sscanf(line12,repmat('%d ',1,numfields));
    end
    line_count = line_count + 1;
    
    % For each field, read the unit and fill value into the data structure
    field_array = cell(1,numfields+1); field_array{1} = 'UTC';
    for f = 1:numfields
        line = fgetl(fid); line_count = line_count + 1;
        comma_pos = strfind(line,',');
        if isempty(comma_pos); %Some files don't comma separate
            comma_pos = strfind(line,' ');
            field = line(1:(comma_pos-1));
            unit = line((comma_pos+1):end);
        else
            field = line(1:(comma_pos-1));
            unit = line((comma_pos+2):end);
        end
        field_array{f+1} = field;
        
        % Sanitize the field name, characters that don't belong in variable
        % names will be removed. Also, prepend 'f_' to any field names 
        field = regexprep(field,'\W','');
        if ~isempty(regexp(field(1),'[^a-zA-Z]')); field = ['f_',field]; end % regexp(field,'[^a-zA-Z]','ONCE') DOES NOT WORK: this line needs to test if the first character in 'field' is not a letter.
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
        if ~isempty(regexp(field(1),'[^a-zA-Z]')); field = ['f_',field]; end % regexp(field,'[^a-zA-Z]','ONCE') DOES NOT WORK: this line needs to test if the first character in 'field' is not a letter.
        header{f} = field;
    end
    
    % Close the fid; we will use dlmread now and no longer need this file
    % open
    fclose(fid);
    
    if DEBUG_LEVEL > 0; fprintf('   Loading the data table\n'); end
    % Read in the full table
    DataTable = dlmread(fullfile(icart_dir, merge_files(a).name),delim,N(1),0);
    
    % Double check that the table size matches the header row
    if size(DataTable,2) ~= numel(header);
        error('read_merge:data_import','Table imported and header row do not have the same number of columns');
    end
    
    % Append the table data to the Merge structure
    for r = 1:numfields+1
        new_field = header{r}; old_field = field_array{r};
        if field_bool(r) == 1;
            field = old_field;
        elseif field_bool(r) == 0;
            field = new_field;
        elseif field_bool(r) == -1;
            if strcmp(new_field,old_field); 
                field_bool(r) = 1;
                field = new_field;
            else
                question = sprintf('Header field %s does not match expected field %s. Which should be used?',new_field, old_field);
                choice = questdlg(question,'Field Mismatch','Header','Existing','Header for all fields','Existing for all fields','Header');
                switch choice
                    case 'Header'
                        field = new_field; field_bool(r) = 0;
                    case 'Existing'
                        field = old_field; field_bool(r) = 1;
                    case 'Header for all fields'
                        field = new_field; field_bool(:) = 0;
                    case 'Existing for all fields'
                        field = old_field; field_bool(:) = 1;
                end
            end
        end
        if DEBUG_LEVEL > 1; fprintf('\tField: %s', field); end
        data_column = DataTable(:,r)';
        eval(sprintf('Merge.Data.%s.Values = data_column;',field));
        if DEBUG_LEVEL > 1; fprintf(' Done.\n'); end
    end
    
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
