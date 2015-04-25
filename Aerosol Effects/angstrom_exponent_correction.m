function [ ext_440nm ] = angstrom_exponent_correction( Merge, campaign_name, meas_wavelength )
%angstrom_exponent_correction Calculate extinction at 440 nm
%   Since aerosols have different optical depths at different wavelengths,
%   it may be important to correct the extinction from the value measured
%   at 532 nm during DISCOVER to the value appropriate for 440 nm (i.e. the
%   band NO2 is typically measured in by OMI).
%
%   This correction is complicated by the fact that scattering and
%   absorption have different behaviors relative to wavelength, so the
%   LARGE group reports their angstrom exponents separately.  Further,
%   these measurements are usually made for dry aerosols (opposed to
%   ambient aerosols) which is to say the aerosols are dried out to remove
%   any water content. Since the satellite is obviously being affected by
%   ambient aerosols - water and all - this is suboptimal.
%
%   The correction follows the formula:
%
%       tau_lambda / tau_0 = (lambda / lambda_0)^(-alpha)
%
%   where alpha is the angstrom exponent.  This function assumes that this
%   relation holds if tau is replaced with k_abs or k_scat, i.e. that the
%   ratio is between an absorption or scattering coefficient, rather than
%   an optical depth.
% 
%   Inputs to this function: 1) a Merge structure, 2) the campaign name (to
%   get field names) and an optional 3) the measurement wavelength.  This
%   is assumed to be 532 nm if unspecified.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VERIFICATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,3);
if ~isstruct(Merge) || ~all(ismember({'metadata','Data'},fieldnames(Merge)))
    % Check that the Merge input is a Merge structure output from
    % read_icart_data.  It should have the fields metadata and Data.
    E.badinput('Merge must be a merge-type structure with fields metadata and Data');
end

if ~ischar(campaign_name)
    E.badinput('campaign_name must be a string');
end

if nargin < 3
    meas_wavelength = 532;
elseif ~(isscalar(meas_wavelength) && isnumeric(meas_wavelength))
    E.badinput('meas_wavelength must be a scalar value');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Wavelength correction %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import data
Names = merge_field_names(campaign_name);
scat_532nm = remove_merge_fills(Merge,Names.aerosol_scattering);
ext_532nm = remove_merge_fills(Merge,Names.aerosol_extinction);
abs_532nm = ext_532nm - scat_532nm;

% Angstrom exponents
ae_abs = remove_merge_fills(Merge,Names.abs_angstr_exp);
ae_scat = remove_merge_fills(Merge,Names.scat_angstr_exp);

abs_440nm = abs_532nm .* (440 / meas_wavelength).^(-ae_abs);
scat_440nm = scat_532nm .* (440 / meas_wavelength).^(-ae_scat);
    
ext_440nm = abs_440nm + scat_440nm;

end

