function [ profile_struct ] = categorize_aerosol_profile(  )
%categorize_aerosol_profile Classify the relation of aerosol profile shape to NO2
%   Bousserez 2014 (ACPD, in preparation) modeled the effect of coincident
%   aerosol and NO2 layers on the AMF for satellite retireved NO2 columns.
%   His models suggest that when part of the NO2 profile lies above the
%   aerosol layer, the AMF should be greater (suggesting and enhancement of
%   the observed SCD by the underlying aerosol layer).  Conversely, when
%   part of the aerosol layer lies above the NO2 layer, he expects that the
%   AMF should be decreased because the aerosol "shields" the lower NO2.
%
%   To examine this effect with in-situ aircraft data, the first step is to
%   establish which profiles fall into which of these categories.  This
%   will do so by comparing the altitude z at which integrating from the
%   surface to z includes 90% of the total column. The results will be
%   placed in 3 categories: (a) those where z_aer and z_NO2 are within dz
%   of each other; (b) those where z_aer > z_NO2 + dz; (c) those where
%   z_NO2 > z_aer + dz.  Profiles will also be further divided into those
%   that have a small maximum aerosol extinction and those with large
%   maximum extinction.  This is adjustable as the "mag_crit" variable.
%
%   All information is returned in a structure: three of the fields contain
%   the criteria used in separating the columns (keeping a record of what
%   criteria led to the given separation) and the last six are the six
%   classes of NO2/aerosol column interation.

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

campaign_name = 'seac4rs';

[Names, dates, directory, range_file] = merge_field_names(campaign_name);

start_date = '';
end_date = '';

if isempty(start_date) || isempty(end_date)
    fprintf('Setting start and end dates based on the campaign name.\n');
    start_date = dates{1};
    end_date = dates{2};
end

%crit_frac is the fraction of the integrated column that we look for to
%compare heights
crit_frac = 0.9;
%mag_crit is the criterion to separate aerosol profiles based on the
%magnitude of the profile.  Set mag_crit_type to 'max' to use the maximum
%extinction in Mm^(-1) or 'int' to use the integrated extinction, i.e.
%column optical depth
mag_crit_type = 'int';
mag_crit = 0.2;
% dz is the minimum separation that the critical altitudes must have to be
% counted as different.
dz = 0.1;

aerosol_field = '';
no2_field = '';
alt_field = '';
alt_conversion = 1; % The navigation data provided with the NO2 merge is already in km usually
profnum_field = '';

merge_directory = '';
merge_file_pattern = '*_%s_%s_%s.mat';


DEBUG_LEVEL = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VARIABLE PREP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
E = JLLErrors;

dates = datenum(start_date):datenum(end_date);
first_warning_alt_unit = 1;

profile_struct = struct('Column_Fraction_Criterion',crit_frac,'Magnitude_Criterion_Type',mag_crit_type,'Magnitude_Critical_Value',mag_crit,'dz',dz,...
    'CoincidentLow',[],'CoincidentLowAlt',[],'CoincidentLowMag',[],'CoincidentHigh',[],'CoincidentHighAlt',[],'CoincidentHighMag',[],'AerosolAboveLow',[],...
    'AerosolAboveLowAlt',[],'AerosolAboveLowMag',[],'AerosolAboveHigh',[],'AerosolAboveHighAlt',[],'AerosolAboveHighMag',[],'NO2AboveLow',[],...
    'NO2AboveLowAlt',[],'NO2AboveLowMag',[],'NO2AboveHigh',[],'NO2AboveHighAlt',[],'NO2AboveHighMag',[]);

% Check that mag_crit_type is one of the two allowed values
if ~any(strcmpi(mag_crit_type,{'max','int'}))
    error(E.badinput('''mag_crit_type'' must be ''max'' or ''int'''));
end

% If the field names are left empty, use the field names returned from
% merge_field_names
if isempty(merge_directory)
    fprintf('Setting merge file directory based on the campaign name.\n');
    merge_directory = directory;
end

if isempty(aerosol_field)
    fprintf('Setting aerosol extinction field based on the campaign name.\n');
    aerosol_field = Names.aerosol_extinction;
end
if isempty(no2_field)
    fprintf('Setting NO2 field based on the campaign name.\n');
    no2_field = Names.no2_lif;
end
if isempty(alt_field)
    fprintf('Setting altitude field based on the campaign name.\n');
    alt_field = Names.gps_alt;
end
if isempty(profnum_field)
    fprintf('Setting profile number field based on the campaign name.\n');
    profnum_field = Names.profile_numbers;
end


%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

% If the profile numbers field is empty, load the range file.  Make a list
% of the dates in the range file to let us match up the date currently
% being worked on with the correct index of the Ranges structure.
if isempty(profnum_field)
    load(range_file); % loads the Ranges variable, which is a data structure
    range_dates = cellstr(datestr({Ranges.Date},29));
end

for d=1:numel(dates)
    % For each day, try to load the aerosol and NO2 merge files.  If one is
    % not available, skip this day
    
    S_NO2 = wildcard_load(merge_directory, merge_file_pattern, dates(d));
    if isempty(S_NO2);
        if DEBUG_LEVEL > 0; fprintf('No NO2 Merge for %s\n',datestr(dates(d))); end
        continue
    else
        Merge = S_NO2.Merge;
    end
    
    % Load the day's aerosol, NO2, and altitude data
    [ext,~,aer_alt] = remove_merge_fills( Merge, aerosol_field, 'alt', alt_field );
    aer_alt = aer_alt * alt_conversion;
    if ~strcmp('km',Merge.Data.(alt_field).Unit) && first_warning_alt_unit;
        warning('Altitude field for merge may not be in km.  Ensure correct conversion factor is being applied');
        first_warning_alt_unit = 0;
    end
    
    [no2,utc,no2_alt] = remove_merge_fills( Merge, no2_field, 'alt', alt_field );
    no2_alt = no2_alt * alt_conversion;
    no2(no2<0) = NaN; % A negative concentration is clearly non-physical
    
    % Indentify each profile present in a given day and loop through each.
    
    if ~isempty(profnum_field)
        if DEBUG_LEVEL > 0; fprintf('Using profile numbers\n'); end
        % If the profile numbers field name is not an empty string, that
        % means that this campaign does indeed have profile numbers.  So
        % we'll use those to determine profiles.
        profnums = Merge.Data.(profnum_field).Values;
        profile_id = unique(profnums);
        
        % Because ranges (the next option) are defined as an n x 2 matrix,
        % we need our list of profile numbers to be a column so that
        % profile_id(a,:) returns a single number.
        if ~iscolumn(profile_id)
            profile_id = profile_id';
        end
        
        % Profile numbers of 0 indicate that the aircraft was not conducting a
        % spiral at that time and so should be removed.
        profile_id(profile_id<1) = [];
        
        % If there were no profiles, let the user know and skip the loop.
        if isempty(profile_id)
            if DEBUG_LEVEL > 0; fprintf('No profiles detected on %s. Skipping.\n',datestr(dates(d))); end
            continue
        end
    else
        % If the profile numbers field name IS an empty string, assume that
        % we need to use UTC ranges instead.  This will depend upon a range
        % file being available. Ranges are to be given as an n x 2 matrix 
        if DEBUG_LEVEL > 0; fprintf('Using UTC ranges\n'); end
        
        xx = find(strcmp(datestr(dates(d),29),range_dates));
        if isempty(xx);
            error('run_spiral:ranges','No UTC ranges found for %s',curr_date);
        end
        profile_id = Ranges(xx).Ranges;    
    end
    
    for p=1:size(profile_id,1)
        % Find the data for each profile from both merges and bin it to
        % 0.25 km bins.  If there is a valid profile number field, find
        % the profile by number, otherwise use the UTC ranges.
        if ~isempty(profnum_field)
            xx = profnums == profile_id(p,:);
        else
            xx = utc >= profile_id(p,1) & utc <= profile_id(p,2);
        end
        
        prof_no2 = no2(xx); prof_no2_alt = no2_alt(xx);
        prof_aer = ext(xx); prof_aer_alt = aer_alt(xx);
        
        [no2_bins, no2_bin_alt] = bin_vertical_profile(prof_no2_alt, prof_no2, 0.25);
        [aer_bins, aer_bin_alt] = bin_vertical_profile(prof_aer_alt, prof_aer, 0.25);
        
        % If only 1 or no bins are not NaNs then we should also skip this
        % profile; there won't be enough information to carry out the
        % analysis
        if allbutone(isnan(no2_bins)) || allbutone(isnan(aer_bins));
            if DEBUG_LEVEL > 1; fprintf('\tNO2 or Aerosol bins for profile #%d on %s have 0 or 1 non-NaN bins\n',profile_id(p,1),datestr(dates(d))); end
            continue 
        end
        
        % Otherwise, trim the leading and trailing NaNs and linearly
        % interpolate any middle NaNs
        [no2_bin_alt, no2_bins] = fill_nans(no2_bin_alt, no2_bins);
        [aer_bin_alt, aer_bins] = fill_nans(aer_bin_alt, aer_bins);
        
        % Now classify the profile: we'll do the full integration for each
        % profile, then integrate over progressively fewer bins until the
        % partial integrated column is < 90% of the total.  By assuming
        % that the integrated column varies quadratically between bin
        % edges we can then solve for the approximate altitude at which the
        % 90% threshold is crossed for both.
        
        % First, get the total integrated profile:
        no2_total = trapz(no2_bin_alt, no2_bins);
        aer_total = trapz(aer_bin_alt, aer_bins);
        
        % Then progressively lower the top altitude until int_z / int_tot <
        % 0.9.  Should the index (x_zi) be 1, that means that the 90% level
        % must be in the first bin, so exit the while loop since trapz()
        % will get confused with the scalar second input.
        no2_partial = no2_total;
        no2_zi = numel(no2_bins);
        while no2_partial / no2_total > crit_frac
            no2_zi = no2_zi - 1;
            if no2_zi == 1; break; end
            no2_partial = trapz(no2_bin_alt(1:no2_zi), no2_bins(1:no2_zi));
        end
        
        aer_partial = aer_total;
        aer_zi = numel(aer_bins);
        while aer_partial / aer_total > crit_frac;
            aer_zi = aer_zi - 1;
            if aer_zi == 1; break; end
            aer_partial = trapz(aer_bin_alt(1:aer_zi), aer_bins(1:aer_zi));
        end
        
        % Now it's likely that the actual point at which the criteria (90%
        % of the column) is met is somewhere in between two bin centers.
        % Since we use the trapezoidal approximation, we can solve for the
        % approximate point where the criteria is met.  The final equation
        % is quadratic, so we need a messy quadratic formula, unless the
        % slope happens to be 0.
        
        % Need the slope between the bins below and above the criteria
        % point
        m_no2 = (no2_bins(no2_zi+1) - no2_bins(no2_zi)) / (no2_bin_alt(no2_zi+1) - no2_bin_alt(no2_zi));
        m_aer = (aer_bins(aer_zi+1) - aer_bins(aer_zi)) / (aer_bin_alt(aer_zi+1) - aer_bin_alt(aer_zi));
        
        if m_no2 == 0
            z90_no2 = (crit_frac * no2_total - no2_partial + no2_bins(no2_zi) * no2_bin_alt(no2_zi))/no2_bins(no2_zi);
        else
            % Define a, b, and c for the quadractic formulae
            a_no2 = 0.5 * m_no2;
            b_no2 = -m_no2 * no2_bin_alt(no2_zi) + no2_bins(no2_zi);
            c_no2 = 0.5 * m_no2 * no2_bin_alt(no2_zi)^2 - no2_bins(no2_zi)*no2_bin_alt(no2_zi) + no2_partial - crit_frac * no2_total;
            
            % Solve for the heights at which 90% of the column is below. Since
            % there are two roots, find the one that is between the bins above
            % and below 90% of the column
            z90_no2_1 = (-b_no2 + sqrt(b_no2^2 - 4*a_no2*c_no2)) / (2 * a_no2);
            z90_no2_2 = (-b_no2 - sqrt(b_no2^2 - 4*a_no2*c_no2)) / (2 * a_no2);
            test = [z90_no2_1, z90_no2_2];
            xx = test > no2_bin_alt(no2_zi) & test < no2_bin_alt(no2_zi+1);
            if imag(z90_no2_1) || imag(z90_no2_2);
                error(E.callError('imaginary','NO2 root is imaginary'));
            elseif all(~xx)
                warning('No solution to NO2 profile. Profile will not be categorized.');
                continue
            elseif all(xx);
                warning('Both solutions for NO2 profile #%d (%0.2f, %0.2f) fall in expected range. Profile will not be categorized.',profile_id(p,1),z90_no2,z90_no2_neg);
                continue;
            else
                z90_no2 = test(xx);
            end
        end
        
        if m_aer == 0
            z90_aer = (crit_frac * no2_total - no2_partial + no2_bins(no2_zi) * no2_bin_alt(no2_zi))/no2_bins(no2_zi);
        else
            % Define a, b, and c for the quadractic formulae
            a_aer = 0.5 * m_aer;
            b_aer = -m_aer * aer_bin_alt(aer_zi) + aer_bins(aer_zi);
            c_aer = 0.5 * m_aer * aer_bin_alt(aer_zi)^2 - aer_bins(aer_zi) * aer_bin_alt(aer_zi)  + aer_partial - crit_frac * aer_total;
            
            % Solve for the heights at which 90% of the column is below. Since
            % there are two roots, find the one that is between the bins above
            % and below 90% of the column
            z90_aer_1 = (-b_aer + sqrt(b_aer^2 - 4*a_aer*c_aer)) / (2 * a_aer);
            z90_aer_2 = (-b_aer - sqrt(b_aer^2 - 4*a_aer*c_aer)) / (2 * a_aer);
            test = [z90_aer_1, z90_aer_2];
            xx = test > aer_bin_alt(aer_zi) & test < aer_bin_alt(aer_zi+1);
            if imag(z90_aer_1) || imag(z90_aer_2);
                error(E.callError('imaginary','Aerosol root is imaginary'));
            elseif all(~xx)
                warning('No solution to aerosol profile. Profile will not be categorized.');
                continue
            elseif all(xx) > 0;
                warning('Both solutions for aerosol profile #%d (%0.2f, %0.2f) fall in expected range. Profile will not be categorized.',profile_id(p,1),z90_aer,z90_aer_neg);
                continue;
            else
                z90_aer = test(xx);
            end
        end
        
        % Set which value to use for assigning the profile to low or high
        % amounts of aerosol
        switch mag_crit_type
            case 'int'
                % In this case, integrate the aerosol profiles taking into
                % account that the units used for extinction are usually
                % Mm^-1 (and yes, that is inverse megameters) and the
                % altitude bins are in km.
                alt_meters = aer_bin_alt * 1000; % convert km to m
                aers_meters = aer_bins * 1e-6; % convert Mm^-1 to m^-1
                
                % Now integrate to get total optical depth of the aerosol
                % column
                aer_opt_depth = trapz(alt_meters,aers_meters);
                
                % Use this as the criterion value for dividing into low and
                % high aerosol loading cases
                mag_value = aer_opt_depth;
            case 'max'
                % In this case we're using the maximum extinction value for
                % a single bin.
                mag_value = max(aer_bins);
        end
        
        % Assign the profile to a specific category
        if abs(z90_aer - z90_no2) < dz && mag_value < mag_crit;
            profile_struct.CoincidentLow = [profile_struct.CoincidentLow; profile_id(p,:)];
            profile_struct.CoincidentLowAlt = [profile_struct.CoincidentLowAlt; profile_id(p,:), z90_no2, z90_aer];
            profile_struct.CoincidentLowMag = [profile_struct.CoincidentLowMag; profile_id(p,:), mag_value];
        elseif abs(z90_aer - z90_no2) < dz && mag_value >= mag_crit;
            profile_struct.CoincidentHigh = [profile_struct.CoincidentHigh; profile_id(p,:)];
            profile_struct.CoincidentHighAlt = [profile_struct.CoincidentHighAlt; profile_id(p,:), z90_no2, z90_aer];
            profile_struct.CoincidentHighMag = [profile_struct.CoincidentHighMag; profile_id(p,:), mag_value];
        elseif z90_aer > z90_no2 + dz && mag_value < mag_crit;
            profile_struct.AerosolAboveLow = [profile_struct.AerosolAboveLow; profile_id(p,:)];
            profile_struct.AerosolAboveLowAlt = [profile_struct.AerosolAboveLowAlt; profile_id(p,:), z90_no2, z90_aer];
            profile_struct.AerosolAboveLowMag = [profile_struct.AerosolAboveLowMag; profile_id(p,:), mag_value];
        elseif z90_aer > z90_no2 + dz && mag_value >= mag_crit;
            profile_struct.AerosolAboveHigh = [profile_struct.AerosolAboveHigh; profile_id(p,:)];
            profile_struct.AerosolAboveHighAlt = [profile_struct.AerosolAboveHighAlt; profile_id(p,:), z90_no2, z90_aer];
            profile_struct.AerosolAboveHighMag = [profile_struct.AerosolAboveHighMag; profile_id(p,:), mag_value];
        elseif z90_no2 > z90_aer + dz && mag_value < mag_crit;
            profile_struct.NO2AboveLow = [profile_struct.NO2AboveLow; profile_id(p,:)];
            profile_struct.NO2AboveLowAlt = [profile_struct.NO2AboveLowAlt; profile_id(p,:), z90_no2, z90_aer];
            profile_struct.NO2AboveLowMag = [profile_struct.NO2AboveLowMag; profile_id(p,:), mag_value];
        elseif z90_no2 > z90_aer + dz && mag_value >= mag_crit;
            profile_struct.NO2AboveHigh = [profile_struct.NO2AboveHigh; profile_id(p,:)];
            profile_struct.NO2AboveHighAlt = [profile_struct.NO2AboveHighAlt; profile_id(p,:), z90_no2, z90_aer];
            profile_struct.NO2AboveHighMag = [profile_struct.NO2AboveHighMag; profile_id(p,:), mag_value];
        end
    end
end

% Output - if there is no output variable, put a variable in the base
% workspace with the name <CAMPAIGN_NAME>_AerCat
if nargout < 1
    % Make the variable name start with the campaign name but in upper case
    % with all non-letter or number characters removed.
    varname = sprintf('%s_AerCat',upper(regexprep(campaign_name,'\W','')));
    eval(sprintf('%s = profile_struct;',varname));
    putvar(sprintf('%s',varname));
end


end

