function [ prof_table, crit_table ] = tabulate_profile_categories( varargin )
%tabulate_profile_categories Returns profiles sorted by date and number
%   To compare the changes in categorization with different criteria, this
%   will take any number of structures returns from
%   categorize_aerosol_profile in and sort the profiles by date and number,
%   appending the categorization given in each structure.  It will return a
%   table with all of this information. The second table will contain a
%   summary of the criteria set for each experiment.
%
%   If you pass a structure of user-defined categorizations, make it the
%   last one and make the final argument a scalar value > 0.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,Inf);

% Figure out if there is a user-defined categorization structure - it'll
% need handled specially
if ~isstruct(varargin{end}) && isscalar(varargin{end}) && varargin{end} > 0
    user_struct = varargin{end-1};
    varargin = varargin(1:end-2);
    user_struct_bool = true;
else 
    user_struct_bool = false;
end

% Check that only structures were passed
for a=1:numel(varargin)
    if ~isstruct(varargin{a})
        error(E.badinput('sort_profile_categories only accepts structures as inputs'));
    end
end


% The "Info" structure will contain important stats about the first
% structure passed, to check against following structure to ensure that
% they are consistent.

fns = fieldnames(varargin{1});
Info.fieldnames = fns;
data_field_bool = true(size(fns));
Prof_Struct = varargin{1};
for a=1:numel(fns)
    if isscalar(Prof_Struct.(fns{a})) || ischar(Prof_Struct.(fns{a}))
        data_field_bool(a) = false;
    end
end
Info.data_fields = fns(data_field_bool);
Info.info_fields = fns(~data_field_bool);
Info.profnum_fields = {'CoincidentLow','CoincidentHigh','AerosolAboveLow','AerosolAboveHigh','NO2AboveLow','NO2AboveHigh'};
Info.date_fields = strcat(Info.profnum_fields,'Dates');

if user_struct_bool
    % A user defined structure from categorize_aerosol_profiles will ONLY
    % have profnum and date fields, so we just need to separate those out
    user_fns = fieldnames(user_struct);
    xx = ~iscellcontents(regexpi(user_fns,'Dates'),'isempty');
    UserInfo.profnum_fields = user_fns(~xx);
    UserInfo.date_fields = user_fns(xx);
end

for a=2:numel(varargin)
    next_fns = fieldnames(varargin{a});
    if numel(next_fns) ~= numel(Info.fieldnames) || ~all(ismember(next_fns, Info.fieldnames))
        E.badinput('All input structures must have the same fields');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Get all the profile numbers and dates from each structure, identify
% unique profiles, and sort them.
Out.all_profnums = [];
Out.all_dates = {};
for a=1:numel(varargin)
    curr_struct = varargin{a};
    for f=1:numel(Info.profnum_fields)
        if ~isempty(curr_struct.(Info.profnum_fields{f}))
            Out.all_profnums = cat(1,Out.all_profnums,curr_struct.(Info.profnum_fields{f})(:,1));
            Out.all_dates = cat(1,Out.all_dates,curr_struct.(Info.date_fields{f}));
        end
    end
end
if user_struct_bool
    for f=1:numel(UserInfo.profnum_fields)
        if ~isempty(user_struct.(UserInfo.profnum_fields{f}))
            Out.all_profnums = cat(1,Out.all_profnums,user_struct.(UserInfo.profnum_fields{f})(:,1));
            Out.all_dates = cat(1,Out.all_dates,user_struct.(UserInfo.date_fields{f}));
        end
    end
end

[Out.all_profnums, sort_perm] = sort(Out.all_profnums);
Out.all_dates = Out.all_dates(sort_perm);
all_datenums = datenum(Out.all_dates);

keep_prof = true(size(Out.all_profnums));

for b=2:numel(Out.all_profnums)
    if Out.all_profnums(b) == Out.all_profnums(b-1) && all_datenums(b) == all_datenums(b-1)
        keep_prof(b) = false;
    end
end

Out.all_profnums = Out.all_profnums(keep_prof);
Out.all_dates = Out.all_dates(keep_prof);
all_datenums = all_datenums(keep_prof);


% Now loop through each structure and make the table. First column is the
% dates, second is prof nums, third and on will be the categorization of
% each profile.

if user_struct_bool
    Out.UserDefined = cell(size(Out.all_profnums));
    for f=1:numel(UserInfo.profnum_fields)
        xx = ismember(Out.all_profnums,user_struct.(UserInfo.profnum_fields{f})) & ismember(all_datenums,datenum(user_struct.(UserInfo.date_fields{f})));
        Out.UserDefined(xx) = UserInfo.profnum_fields(f);
    end
end

for a=1:numel(varargin)
    curr_struct = varargin{a};
    curr_table_col = sprintf('InputStct%02d',a);
    Out.(curr_table_col) = cell(size(Out.all_profnums));
    for f=1:numel(Info.profnum_fields)
        xx = ismember(Out.all_profnums,curr_struct.(Info.profnum_fields{f})) & ismember(all_datenums,datenum(curr_struct.(Info.date_fields{f})));
        Out.(curr_table_col)(xx) = Info.profnum_fields(f);
    end
    
    Criteria.(curr_table_col) = cell(numel(Info.info_fields),1);
    for f=1:numel(Info.info_fields)
        Criteria.(curr_table_col){f} = curr_struct.(Info.info_fields{f});
    end
end

% Convert both structures to tables, and in the case of the Criteria table,
% assign the row variable names.

prof_table = struct2table(Out);
crit_table = struct2table(Criteria,'RowNames',Info.info_fields);


end

