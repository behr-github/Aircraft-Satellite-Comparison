function matches = verify_sat_vs_pandora(Merge, Data, varargin)
%VERIFY_SAT_VS_PANDORA Verify satellite data against Pandora columns
%   MATCHES = VERIFY_SAT_VS_PANDORA( MERGE, DATA ) Matches Pandora total
%   NO2 columns in MERGE against OMI SP or BEHR VCDs (with adding SP
%   stratospheric NO2). MERGE must be a Merge structure generated by
%   reading Pandora ICARTT files from one of the DISCOVER campaigns by
%   READ_ICARTT_FILES(), or a structure with the same data hierarchy. DATA
%   is a BEHR Data structure.
%
%   Parameters:
%       'time_range' - time in hours to either side of the OMI overpass
%       time for matched pixels to average Pandora data for comparison to
%       OMI column data.


% In contrast to verify_sat_vs_aircraft_vcds(), because Pandora instruments
% sit in one place and take measurements continuously, we need to first
% match up which satellite pixels the Pandora instrument falls in. Once
% we've done that, it'll be easier to match the Pandora times since we'll
% know the exact time of each OMI pixel.

matches = match_pixels_to_pandora(Merge, Data);
if isempty(matches)
    matches = format_output(matches);
    return
end

matches = match_pandora_to_overpass(Merge, matches, varargin{:});
matches = format_output(matches);

end

%%%%%%%%%%%%%%%%%%%%%
% MAIN SUBFUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%

function matches = match_pixels_to_pandora(Merge, Data)
% In order to match OMI pixels to Pandora sites, first we need to get the
% Pandora location. 

[pandora_lon, pandora_lat] = get_pandora_lat_lon(Merge);

% Find all pixels that contain this lat/lon in each orbit. Unlike the
% profiles, we can't check for valid measurements yet. That will have to
% come after we match the times.
matches = [];
for i_orbit = 1:numel(Data)
    xx_pix = pixels_overlap_profile(Data(i_orbit).FoV75CornerLongitude, Data(i_orbit).FoV75CornerLatitude, pandora_lon, pandora_lat);
    
    if sum(xx_pix(:)) == 0
        continue
    end
    
    % If we found a or pixels, reject invalid pixels, then check if any of
    % the pixel matched are valid.
    Data(i_orbit).Areaweight = ones(size(Data(i_orbit).Longitude));
    xx_good_pix = filter_bad_pixels(Data(i_orbit));
    
    if sum(xx_pix(:) & xx_good_pix(:)) == 0
        continue
    end
    
    if isempty(matches)
        matches = get_valid_matched_data(Data(i_orbit), xx_good_pix & xx_pix, pandora_lon, pandora_lat);
    else
        matches = cat_fields(matches, get_valid_matched_data(Data(i_orbit), xx_good_pix & xx_pix, pandora_lon, pandora_lat));
    end
end

end


function matches = match_pandora_to_overpass(Merge, matches, varargin)
E = JLLErrors;
p = inputParser;
p.addParameter('time_range',1);
p.KeepUnmatched = true;

p.parse(varargin{:});
pout = p.Results;

time_range = pout.time_range/24; % convert from hours to days so that it can be used directly with datenums

pandora_dnums = pandora_gmt_to_datenum(remove_merge_fills(Merge, 'DateGMT'));
pandora_no2 = remove_merge_fills(Merge, pandora_no2_fieldname(Merge));

xx_valid = filter_good_pandora_obs(Merge);

% This should ensure that if there's no Pandora data near enough in time to
% the OMI overpass that there will be a NaN as a fill value.
matches.pandora_no2 = nan(size(matches.omi_time));

for i_match = 1:numel(matches.omi_time)
    xx_time = abs(pandora_dnums - matches.omi_time(i_match)) < time_range;
    matches.pandora_no2(i_match) = nanmean(pandora_no2(xx_valid & xx_time));
    matches.details(i_match).all_pandora_no2 = pandora_no2(xx_valid & xx_time);
    matches.details(i_match).all_pandora_dnums = pandora_dnums(xx_valid & xx_time);
end

end

%%%%%%%%%%%%%%%%%%%%%%%
% HELPER SUBFUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%

function xx_good = filter_bad_pixels(Data)
Data = omi_pixel_reject(Data, 'detailed', struct('cloud_type', 'omi', 'cloud_frac', 0.2, 'row_anom_mode', 'XTrackFlags', 'check_behr_amf', true));
xx_good = Data.Areaweight > 0;
end

function matches = get_valid_matched_data(Data, xx_match, pandora_lon, pandora_lat)

matches.sp_no2 = nanmean(Data.ColumnAmountNO2Trop(xx_match) + Data.ColumnAmountNO2Strat(xx_match));
matches.behr_no2 = nanmean(Data.BEHRColumnAmountNO2Trop(xx_match) + Data.ColumnAmountNO2Strat(xx_match));
time_2d = repmat(Data.Time, 1, size(Data.Longitude,2));
matches.omi_time = omi_time_conv(nanmean(time_2d(xx_match)));
matches.details = subset_behr_data_struct(Data, xx_match);
% Name these "profile_lon" and "profile_lat" to match the naming in the
% aircraft structure.
matches.profile_lon = repmat(pandora_lon, size(matches.sp_no2));
matches.profile_lat = repmat(pandora_lat, size(matches.sp_no2));

end

function matches = format_output(matches)
% Ensure that matches always has the same fields when output. This way if
% we want to concatenate the match structures we can.
req_fields = {'sp_no2', 'behr_no2', 'omi_time', 'details', 'pandora_no2', 'profile_lon', 'profile_lat'};
if isempty(matches)
    matches = make_empty_struct_from_cell(req_fields);
else
    fns_matches = fieldnames(matches);
    missing_fields = req_fields(~ismember(req_fields, fns_matches));
    if numel(missing_fields) > 0
        % If there are missing fields then it is possible that if we try to
        % correlate say sp_no2 and pandora_no2 from multiple calls to this
        % main function that there will be different numbers of
        % observations, which would be bad, because it would throw off
        % which observations are matched. So just error if we're missing
        % fields, don't try to assume anything about what they should be.
        E.callError('missing_field', 'matches is missing the fields ''%s'' at output', strjoin(missing_fields, ''', '''));
    end
    
    % If there are fields in matches but not in req_fields, issue a warning
    xx_new = ~ismember(fns_matches, req_fields);
    if any(xx_new)
        warning('Extra fields in matches (%s) not defined in req_fields. Update req_fields', strjoin(fns_matches(xx_new), ', '));
    end
end
end

function xx_good = filter_good_pandora_obs(Merge)
% Returns indices of NO2 observation that should be cloud free and good
% quality. As recommended in the MD Aldino file, filter for Norm_RMS < 0.05
% and Error < 0.005 (strict). These criteria appear to already be applied
% for the Filtered files, this function ensures that if raw files are used
% for some reason, bad data is not included.
norm_rms = remove_merge_fills(Merge,'Norm_RMS');
no2_error = remove_merge_fills(Merge,'Error');
xx_good = norm_rms < 0.05 & no2_error < 0.005;
end