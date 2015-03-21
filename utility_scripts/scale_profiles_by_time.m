function [ output_args ] = scale_profiles_by_time( prof_times, campaign_name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


campaign_name = lower(regexprep(campaign_name,'[\W_]',''));

switch campaign_name
    case 'discovermd'
        bintimes = [7.5000 9.5000 10.5000 11.5000 12.5000 13.5000 14.5000 15.5000 16.5000];
        binscale = [1.1648 0.9422 1.0607 1.0771 0.8195 1.0000 1.1866 1.1956 1.4429];
    otherwise
        

end

