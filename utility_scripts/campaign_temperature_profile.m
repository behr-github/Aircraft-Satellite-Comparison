function [ temperature_prof, pres_profile ] = campaign_temperature_profile( campaign_name )
%campaign_temperature_profile Returns the median temp profile for a given campaign.
%   Over the course of a typical field campaign, it is unlikely that the
%   free troposphere temperature will vary that significantly, so these
%   profiles were derived from temperature data for entire campaign (see
%   make_composite_profile.m in NO2 Profiles/One-off scripts).

E = JLLErrors;

if ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'md'))
    temperature_prof = [300.4550 299.9400 298.9150 298.9100 299.1000 299.0700 296.6100 295.0600 293.5300 290.4100 290.9900 287.9900 286.7000 285.6500 284.5300 281.0800 280.9000 276.4500 272.8400 268.0562 263.6088 258.6371 253.0005 243.5814 229.3784];
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'ca'))
    temperature_prof = [283.5403 284.4901 283.9849 284.1614 294.1300 296.7400 287.3189 292.4100 288.3163 287.3151 288.9400 286.2400 283.9400 282.3737 279.8776 280.3500 280.7500 276.2100 272.7200 251.4486 251.5334 251.6282 251.7357 251.9153 252.1861];
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'tx'))
    temperature_prof = [296.5900 299.4791 288.6492 287.3828 297.9100 297.4400 293.9000 294.5200 291.6400 288.9400 289.9278 287.5900 285.8800 284.1500 280.3100 281.1402 280.0400 276.4336 272.7800 251.4486 251.5334 251.6282 251.7357 251.9153 252.1861];
elseif ~isempty(regexpi(campaign_name,'seac4rs')) || ~isempty(regexpi(campaign_name,'seacers'));
    temperature_prof = [300.6500 300.2505 297.6500 288.3628 298.1900 297.6500 296.1088 294.9000 293.0738 290.9000 290.3594 288.6500 286.6884 285.1200 280.7724 281.5697 279.3009 275.6400 272.0000 265.1500 262.6500 254.1500 248.9000 237.9000 224.9000];
elseif ~isempty(regexpi(campaign_name,'intex')) && ~isempty(regexpi(campaign_name,'b'))
    temperature_prof = [300.2849 287.8000 289.8130 289.1352 297.9000 297.4000 295.9000 294.7948 292.7500 290.8000 290.1500 288.4824 286.3800 284.9700 280.9400 281.4000 279.0000 274.9000 270.9000 263.4000 261.0000 252.2000 243.5000 232.4000 221.1500];
else
    error(E.badinput('Could not parse the given campaign name'));
end

pres_profile = bin_omisp_pressure(1,1);
if temperature_prof(end) > 230
    warning('Temperature profile does not reach the expected range of temperatures for the tropopause. Use with caution.');
end
    
end

