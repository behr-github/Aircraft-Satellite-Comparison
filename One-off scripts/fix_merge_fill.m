% set_merge_fill
%
%   Some of the first Merge files I imported do not have a fill value set
%   for the UTC field, which messes things up in other scripts.  This will
%   correct that for.

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
merge_file_pattern = '*.mat';
DEBUG_LEVEL = 2;

files = dir(fullfile(merge_dir, merge_file_pattern));
for a=1:numel(files)
    filename = files(a).name;
    if DEBUG_LEVEL > 0; fprintf('Loading %s\n',filename); end
    load(fullfile(merge_dir,filename));
    if isfield(Merge.Data.UTC,'fill');
        if DEBUG_LEVEL > 1; fprintf('\tRemoving "fill"\n'); end
        Merge.Data.UTC = rmfield(Merge.Data.UTC, 'fill');
    end
    if ~isfield(Merge.Data.UTC,'Fill');
        if DEBUG_LEVEL > 1; fprintf('\tAdding "Fill"\n'); end
        Merge.Data.UTC.Fill = NaN;
    end
    if DEBUG_LEVEL > 0; fprintf('\tSaving %s\n',filename); end
    save(fullfile(merge_dir,filename),'DataTable','Merge','header');
    
end