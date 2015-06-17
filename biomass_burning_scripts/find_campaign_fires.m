function [ varargout ] = find_campaign_fires( Merge )
%calc_campaign_background Calculate day-by-day background concentrations
%   Detailed explanation goes here

E = JLLErrors;
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% Which campaign to calculate the background for
campaign_name = 'seac4rs';

% Data fields of interest to be collected from each fire plume. These
% should correspond to valid fields in the Names structure returned from
% merge_field_names, e.g. 'no2_lif', 'co', 'aerosol_extinction'
%data_fields = {'no2_lif', 'aerosol_extinction', 'pressure', 'longitude', 'latitude'};
data_fields = {};

% Use a UTC ranges file? Otherwise will just look for enhanced CO.
use_utc_ranges = false;

% The averaging period for enhancement data. Higher values will
% probably lead to more contiguous plumes, but can miss shorter ones.
% Alvarado 2010 (ACP p. 9739) used 60 s average CO data. Note that
% unaveraged data is used for the correlation test.
n_sec_avg = 15;

% If not using range files, this will be how many averaging periods
% separate fire blocks can be apart and be considered 1 plume. E.g. if
% n_sep is 1, and periods 1-4, 6-7, and 11-15 are considered fire plumes,
% then 1-4 will be joined with 6-7 (and all data from minute 1-7 will be
% used) but 11-15 will be a separate plume.
n_sep = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[Names, ~, merge_dir] = merge_field_names(campaign_name);


% Check that the data fields requested are present in the Names structure
data_chk = ismember(data_fields, fieldnames(Names));
if any(~data_chk)
    bad_fields = strjoin(data_fields(~data_chk),', ');
    E.badinput('The fields %s are not defined in the Names structure',bad_fields);
else
    data_chk = false(size(data_fields));
    for a=1:numel(data_fields)
        if isempty(Names.(data_fields{a}))
            data_chk(a) = true;
        end
    end
    if any(~data_chk)
        bad_fields = strjoin(data_fields(~data_chk),', ');
        E.badinput('The fields %s are empty for the campaign %s', bad_fields, campaign_name);
    end
end


% If the user wants to use a range file, load it now.  We do this
% separately from the other call to merge_field_names because if there is
% >1 range file, it will prompt the user to select one.
if use_utc_ranges
    [~,~,~,range_file] = merge_field_names(campaign_name);
    RF=load(range_file);
    Ranges=RF.Ranges;
    range_dates = datenum({Ranges.Date}); %we'll need this to find the right set of ranges
end

% This will be the output variable, a structure containing UTC ranges of
% fire plumes for each day
Fires = struct;

if nargin < 1
    % Get all the Merge files for this campaign, unless the user passed a
    % specific Merge
    F = dir(fullfile(merge_dir,'*.mat'));
    nF = numel(F);
else
    nF = 1;
end
% Load each file in sequence.
merge_dates = cell(nF,1);
for a=1:nF
    if nargin < 1
        load(fullfile(merge_dir,F(a).name),'Merge');
    end
    
    %Check the CO, ACN, and HCN data that they are not all fill values (or
    %generall < 0, since that doesn't make sense for a concentration)
    if all(Merge.Data.(Names.co).Values < 0) || all(Merge.Data.(Names.acn).Values < 0) || all(Merge.Data.(Names.hcn).Values < 0)
        if DEBUG_LEVEL > 0; fprintf('Skipping %s - too much missing data\n',Merge.metadata.date); end
        continue
    end
    
    merge_dates{a} = Merge.metadata.date;
    
    utc = remove_merge_fills(Merge, 'UTC');
    
    co_enhancement = calculate_enhancement(Merge, 'co', campaign_name);
    [acn_enhancement, acn_bck_stddev] = calculate_enhancement(Merge, 'acn', campaign_name);
    [hcn_enhancement, hcn_bck_stddev] = calculate_enhancement(Merge, 'hcn', campaign_name);
    
    
    % Calculate enhancement values - these will be used to determine the
    % presence of fires.
    co_enh_min_avg = avg_n_elements(co_enhancement,n_sec_avg,'op','nanmean'); %The initial CO, ACN, and HCN data is second data
    acn_enh_min_avg = avg_n_elements(acn_enhancement,n_sec_avg,'op','nanmean'); 
    hcn_enh_min_avg = avg_n_elements(hcn_enhancement,n_sec_avg,'op','nanmean'); 

    utc_min_start_time = utc(1:n_sec_avg:end);
    utc_min_end_time = utc(n_sec_avg:n_sec_avg:end);
    % Cut down the start and end time vectors to be the same length as the
    % CO enhancement vector
    l = length(co_enh_min_avg);
    utc_min_start_time = utc_min_start_time(1:l);
    utc_min_end_time = utc_min_end_time(1:l);
    
    % All three vectors should now be the same length
    if length(co_enh_min_avg) ~= length(utc_min_start_time) || length(co_enh_min_avg) ~= length(utc_min_end_time) || length(co_enh_min_avg) ~= length(acn_enh_min_avg) || length(co_enh_min_avg) ~= length(hcn_enh_min_avg)
        E.badvar('co_enh_min_avg, acn_enh_min_avg, hcn_enh_min_avg, utc_min_start_time, utc_min_end_time','These are not all the same length')
    end
    
    if isempty(regexp(Merge.Data.(Names.co).Unit,'ppb','ONCE'))
        % Check that the units of CO are in ppb or ppbv
        E.callError('badunit','Expecting CO data to be in ppb - cannot confirm');
    elseif isempty(regexp(Merge.Data.(Names.acn).Unit,'ppt','ONCE'))
        % Check that the units of ACN are in ppt or pptv
        E.callError('badunit','Expecting ACN data to be in ppt - cannot confirm');
    elseif isempty(regexp(Merge.Data.(Names.hcn).Unit,'ppt','ONCE'))
        % Check that the units of HCN are in ppt or pptv
        E.callError('badunit','Expecting HCN data to be in ppt - cannot confirm');
    end
   
    
    if use_utc_ranges
        rr = range_dates == datenum(Merge.metadata.date);
        if sum(rr) == 0
            E.callError('date_invalid','Could not find UTC range data for %s', Merge.metadata.date);
        elseif sum(rr) > 1
            E.callError('date_invalid','Multiple UTC ranges identified for %s', Merge.metadata.date);
        end
        
        is_fire = false(size(Ranges(rr).Ranges,1),1);
        for b=1:numel(is_fire)
            block_start_utc = Ranges(rr).Ranges(b,1);
            block_end_utc = Ranges(rr).Ranges(b,2);
            is_fire(b) = is_block_fire(block_start_utc, block_end_utc);
        end
        
        fire_ranges = Ranges(rr).Ranges(is_fire,1:2);
    else
        co_enh_bool = double(co_enh_min_avg > 20); % This criterion of 20 ppb above background comes from Alvarado (2010)
        % Since Alvarado didn't have a specific number they looked for for
        % ACN/HCN enhancement (just that it was correlated with CO
        % enhancement), let's take (for now) 20 ppt - consistent with the
        % units used, at least.
        acn_hcn_enh_bool = 2*double(acn_enh_min_avg > 20 | hcn_enh_min_avg > 20);
        
        % This will ID separately areas with no enhancement, CO only,
        % ACN/HCN only, and CO and ACN/HCN. This should help deal with
        % cases where anthropogenic and biomass burning plumes are close to
        % each other temporally.
        blocks = findBlock(co_enh_bool + acn_hcn_enh_bool);
        
        % Now we only want times when CO is enhanced, which means that the
        % third column of blocks will be odd (since a value of 1 was given
        % above for CO enhanced and added to 2 if ACN/HCN enhanced). Set
        % times when only ACN/HCN is enhanced to 0 (not a fire)
        blocks(blocks(:,3)==2,3) = 0;
        
        % For each block of CO enhancement, calculate the correlation between
        % CO enhancement and HCN/ACN enhancement.  If r^2 > 0.3, keep this
        % marked as a fire block (c.f. Alvarado et al. 2010, ACP p. 9739)
        for b=1:size(blocks,1)
            % Test for ACN or HCN enhancement correlation with CO
            % enhancement.
            if blocks(b,3) > 0
                block_start_utc = utc_min_start_time(blocks(b,1));
                block_end_utc = utc_min_end_time(blocks(b,2));
                blocks(b,3) = is_block_fire(block_start_utc, block_end_utc);
            end
        end
        % At this point the final column of blocks will be 1 or 0 (fire/not
        % a fire). Now we want to combine them. First we'll just combine
        % adjacent fire/not fire blocks.
        for b=(size(blocks,1)):-1:2
            if blocks(b,3) == blocks(b-1,3)
                blocks(b-1,2) = blocks(b,2);
                blocks(b,:) = [];
            end
        end
        
        % Allow the user to treat blocks separated by n_sep minutes as one
        % continuous
        if n_sep > 0
            for b=(size(blocks,1)-1):-1:2
                % Since we ran findBlock on a logical matrix, there will be
                % alternating blocks of true and false. To make blocks of
                % true contiguous if they're separated by a false block
                % smaller than the criteria
                if (b < size(blocks,1) && blocks(b,3) == blocks(b+1,3)) || (b>1 && blocks(b,3) == blocks(b-1,3))
                    % This means that the blocks are NOT alternating and something went
                    % very wrong
                elseif ~blocks(b,3) && blocks(b,2) - blocks(b,1) < n_sep;
                    % If it is a false block and it is smaller than n_sep, combine the
                    % blocks before and after. Since we're running backwards, we can do
                    % this in place.
                    blocks(b-1,2) = blocks(b+1,2);
                    blocks(b:b+1,:) = [];
                end
            end
        end
        is_fire = logical(blocks(:,3));
        fire_ranges = blocks(is_fire,1:2);
        % Convert back to UTC values
        for b=1:size(fire_ranges,1)
            fire_ranges(b,1) = utc_min_start_time(fire_ranges(b,1));
            fire_ranges(b,2) = utc_min_end_time(fire_ranges(b,2));
        end
    end
    
    merge_date = regexprep(Merge.metadata.date,'\D','');
    fieldname_date = sprintf('Flight_%s',merge_date);
    Fires.(fieldname_date).utc_ranges = fire_ranges;
    for f=1:numel(data_fields)
        Fires = add_fire_data(Merge, Fires, data_fields{f}, fieldname_date);
    end
end

if nargout > 0
    varargout{1} = Fires;
else
    putvar(Fires);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NESTED FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function fire_bool = is_block_fire(block_start, block_end)
        utc_xx = utc >= block_start & utc < block_end;
        block_co = co_enhancement(utc_xx);
        block_acn = acn_enhancement(utc_xx);
        block_hcn = hcn_enhancement(utc_xx);
        
        [~,~,~,FitData] = calc_fit_line(block_co, block_acn, 'regression', 'rma');
        r2_acn = FitData.R2;
        [~,~,~,FitData] = calc_fit_line(block_co, block_hcn, 'regression', 'rma');
        r2_hcn = FitData.R2;
        
        % If neither ACN or HCN are correlated with the CO enhancement,
        % this is not a fire plume (Alvarado 2010)
        fire_bool = r2_acn > 0.3 || r2_hcn > 0.3;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUBFUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function [enhancement, bckgnd_stddev] = calculate_enhancement(Merge, species_name, campaign_name)
    Names = merge_field_names(campaign_name);
    [species_bckgnd, bckgnd_stddev] = calc_day_background(Merge, species_name, campaign_name);
    
    species = remove_merge_fills(Merge, Names.(species_name));   
    enhancement = species - species_bckgnd;
end

function Fires = add_fire_data(Merge, Fires, field_name, fieldname_date)
[data, utc] = remove_merge_fills(Merge, field_name);
utc_ranges = Fires.(fieldname_date).utc_ranges;
s = size(utc_ranges,1);
Fires.(fieldname_date).(field_name) = cell(s,1);

for r=1:s
    xx = utc >= utc_ranges(r,1) & utc < utc_ranges(r,2);
    Fires.(fieldname_date).(field_name){r} = data(xx);
end
end

