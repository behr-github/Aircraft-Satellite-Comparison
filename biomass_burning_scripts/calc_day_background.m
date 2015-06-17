function [species_background, background_std_dev] = calc_day_background(Merge, species_name, campaign_name, bottom_percent)
% For the given species, calculate the background as the x'th quantile
% of the data. By default, x is 0.1, but this can be overridden with the
% optional fourth argument. Returns the background value and the standard
% deviation of all data in that bottom quantile.
%
%   Josh Laughner <joshlaugh5@gmail.com> 16 June 2015

E=JLLErrors;

Names = merge_field_names(campaign_name);

if ~isfield(Names, species_name)
    E.badinput('%s is not a defined species for %s', species_name, campaign_name)
end

if nargin < 4
    bottom_percent = 0.1;
end

species = remove_merge_fills(Merge, Names.(species_name));
species(isnan(species)) = [];

species_background = quantile(species, bottom_percent);
species_subset = species(species < species_background);
background_std_dev = std(species_subset);

end

