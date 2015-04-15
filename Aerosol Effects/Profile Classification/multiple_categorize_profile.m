function [ OutStruct ] = multiple_categorize_profile( no2_in, no2_alt_in, aer_in, aer_alt_in )
%multiple_categorize_profile Run categorize_aerosol_profile 3 times and take the most common result.
%   The general idea of categorizing the relative position of NO2 and
%   aerosol profiles by comparing the altitude below which x% of the column
%   can be found works well, but can sometimes ignore the structure of the
%   profile at lower altitudes.  By checking the categorization for
%   multiple critical fractions and taking the most common result from
%   them, this should overcome that problem.
%
%   Josh Laughner <joshlaugh5@gmail.com> 7 Apr 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%
%%%%% INPUT %%%%%
%%%%%%%%%%%%%%%%%

campaign_name = 'discover-ca';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL CATEGORIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

crit_fracs = [0.5, 0.75, 0.9];
cat_structs = cell(1,numel(crit_fracs));
for a=1:numel(crit_fracs)
    if nargin < 4
        cat_structs{a} = categorize_aerosol_profile('campaign_name',campaign_name,'crit_frac',crit_fracs(a),'DEBUG_LEVEL',1);
        out_fxn = @add_cat;
    else
        cat_structs{a} = categorize_aerosol_profile(no2_in, no2_alt_in, aer_in, aer_alt_in,'crit_frac',crit_fracs(a),'DEBUG_LEVEL',1);
        out_fxn = @add_cat_simple;
    end
end

if nargin < 4
    % This is the case where we are dealing with an entire campaign's worth
    % of data.
    % Since tabulate_profile_categories already nicely arranges everything but
    % profile we'll go ahead and use that.
    [cat_tab, crit_tab] = tabulate_profile_categories(cat_structs{:});
    
    cat_data = table2cell(cat_tab);
    % Column 1 = profile #, Column 2 = date string, Column 3+ = category
else
    % If, on the other hand, we're just looking at a specific profile, then
    % we've already got a cell array of categorizations.
    cat_data = [cell(1,2), cat_structs];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% META-CATEGORIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prep the output structure: it will have a similar form to the usual
% categorize_aerosol_prof output structures, although it will lose some of
% the extra fields beyond profile ID and date.

if nargin < 4
    OutStruct.Column_Fraction_Criterion = crit_tab{'Column_Fraction_Criterion',:};
    OutStruct.Magnitude_Criterion_Type = crit_tab{'Magnitude_Criterion_Type',:};
    OutStruct.Magnitude_Critical_Value = crit_tab{'Magnitude_Critical_Value',:};
    OutStruct.dz = crit_tab{'dz',:};
end

OutStruct.CoincidentLow = [];
OutStruct.CoincidentLowDates = {};
OutStruct.CoincidentHigh = [];
OutStruct.CoincidentHighDates = {};
OutStruct.AerosolAboveLow = [];
OutStruct.AerosolAboveLowDates = {};
OutStruct.AerosolAboveHigh = [];
OutStruct.AerosolAboveHighDates = {};
OutStruct.NO2AboveLow = [];
OutStruct.NO2AboveLowDates = {};
OutStruct.NO2AboveHigh = [];
OutStruct.NO2AboveHighDates = {};

% Rules: if two or three categorizations agree on the location
% (coinc/aer above/no2 above), that categorization is chosen.

for a=1:size(cat_data,1)
    cats = cat_data(a,:);
    % Look for the categories - allow there to be a space in between
    % aerosol & above or no2 & above.
    xx_coinc = ~iscellcontents(regexpi(cats(3:end),'Coincident'),'isempty');
    xx_aer = ~iscellcontents(regexpi(cats(3:end),'Aerosol ?Above'),'isempty');
    xx_no2 = ~iscellcontents(regexpi(cats(3:end),'NO2 ?Above'),'isempty');
    
    if sum(xx_coinc) > 1
        OutStruct = out_fxn(OutStruct, cats, xx_coinc);
    elseif sum(xx_aer) > 1
        OutStruct = out_fxn(OutStruct, cats, xx_aer);
    elseif sum(xx_no2) > 1
        OutStruct = out_fxn(OutStruct, cats, xx_no2);
    end
end

end

function OutStruct = add_cat(OutStruct, cat_line, xx)
    % find the categorization name - we want to read it from the table to
    % capture whether it is "low" or "high" AOD
    E = JLLErrors;
    cats = cat_line(3:end);
    cat_name = unique(cats(xx));
    if numel(cat_name) > 1 % both names should be the same. if not, don't continue.
        E.badvartype(cat_name,'1x1 cell array');
    end
    
    catdate_name = strcat(cat_name,'Dates');
    catprof = cat_line{1};
    catdate = cat_line(2); % we want the actual cell, not its contents, to concatenate
    
    OutStruct.(cat_name{1}) = cat(1,OutStruct.(cat_name{1}),catprof);
    OutStruct.(catdate_name{1}) = cat(1,OutStruct.(catdate_name{1}),catdate);
end

function cat_name = add_cat_simple(~, cat_line, xx)
    E = JLLErrors;
    cats = cat_line(3:end);
    cat_name = unique(cats(xx));
    if numel(cat_name) > 1 % both names should be the same. if not, don't continue.
        E.badvartype(cat_name,'1x1 cell array');
    end
end

