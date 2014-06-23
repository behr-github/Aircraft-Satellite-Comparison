function [ blh, utc ] = findall_no2_bdy_layer_heights(utc_in,no2,alt_no2,ranges_in)
%findall_no2_bdy_layer_heights(utc, no2, ranges): Finds the boundary layer crossing for each range passed.
%   This function finds the NO2 chemical boundary layer in aircraft data by
%   examining the NO2 data during times when the aircraft is constantly
%   ascending or descending.  These times are passed as the variable
%   ranges, which must be an n x 2 matrix where each row is the start/stop
%   time in seconds after midnight UTC for each range.  Several MAT files
%   with these ranges stored can be found at ~/MATLAB/NO2
%   Profiles/Workspaces
%
%   utc_in and no2 are the values from the aircraft merge file.
%
%   This returns two values, the boundary layer heights and median UTC
%   times of each range.

% Remove any ranges that are entirely outside our UTC range
range_start_test = ranges_in' < min(utc_in);
range_end_test = ranges_in' > max(utc_in);

% As of Matlab v. R2013b, "all" acts on matrices (not vectors) by testing
% if every element in a column (first index) is true. This line thus checks
% if both the start and end of the range are before the lowest UTC value or
% after the highest UTC value.
out_of_range = all(range_start_test) | all(range_end_test);

ranges = ranges_in(~out_of_range',:);

% If no ranges are left, return empty matrices
if isempty(ranges);
    blh = [];
    utc = [];
else
    
    s=size(ranges);
    blh = -99*ones(s(1),1);
    utc = -99*ones(s(1),1);
    
    for a=1:s(1)
        startind = find(abs(utc_in - ranges(a,1)) == min(abs(utc_in - ranges(a,1))));
        endind = find(abs(utc_in - ranges(a,2)) == min(abs(utc_in - ranges(a,2))));
        no2seg = no2(startind:endind);
        altseg = alt_no2(startind:endind);
        utc(a) = nanmedian(utc_in(startind:endind));
        
        % Bin the data so that the vertical profiles are less noisy
        [no2seg, altseg] = bin_vertical_profile(altseg,no2seg,0.5);
        
%         % For most measurements, to find the boundary layer, we look for the
%         % altitude at which the [NO2] is 1/e of its maximum value (as well as
%         % requiring that d[NO2]/dz < 0).  However, when 1 e-fold of the max
%         % value is < 40 pptv, this method tends to run into problems because
%         % the average [NO2] above the BL is ~40 pptv. Thus, when this is the
%         % case, we look for the greatest gradient in NO2, with d[NO2]/dz < 0
%         if max(no2seg(:))*exp(-1) < 40;
%             method = 'max';
%         else
%             method = 'exp';
%         end
        method = 'max'; % 23 Jun 2014: We're going to try just looking for the maximum gradient to see if that brings the column densities down
        blh(a) = find_bdy_layer_height(no2seg, altseg,method);
    end
    
    % Clean up the output: if any fill values (-99) are left, remove those
    % entries
    
    fills = blh == -99;
    blh = blh(~fills);
    utc = utc(~fills);
end
end

