function [ cb ] = plot_flightpath( Merge_str, data_point_skip )
%plot_flightpath(Merge_str): Plots the pressure altitude and flight path from an ICART merge structue
%   This function will automatically plot the flight path of an ICART merge
%   file, coloring the line by altitude.  
%
%   The first argument must be a Merge data structure created by the
%   read_merge_data script.
%
%   The optional second argument will reduce the number of points plotted,
%   e.g. (Merge, 4) would only plot every 4th point.
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 May 2014

p = inputParser;
p.addRequired('Merge_str');
p.addOptional('skip',1,@isscalar);

p.parse(Merge_str, data_point_skip)
pout = p.Results;
Merge= pout.Merge_str;
skip = pout.skip;

% Read in the pressure altitude, lat, and lon data:
p_alt = Merge.Data.PRESSURE.Values;
lat = Merge.Data.LATITUDE.Values;
lon = Merge.Data.LONGITUDE.Values;
lon = lon - 360; %Since m_map uses the convention that W is < 0, we must make the conversion.

% Get the current colormap; we'll need it to color the line
cm = colormap;
cml = size(cm,1); 

% Reduce the size of the plotted vectors
tmp = true(size(lat));
tmp(1:skip:end) = 0;

lat(tmp) = [];
lon(tmp) = [];
p_alt(tmp) = [];

% Calculate the lat/lon boundaries with some padding
latbdy = [floor(min(lat))-0.5, ceil(max(lat))+0.5];
lonbdy = [floor(min(lon))-0.5, ceil(max(lon))+0.5];

% Calculate an integer version of pressure altitude that has its lowest
% value at 1 and its highest as the first dimension size of the colormap.
% This is necessary so that these can serve as indices to the colormap.
p_alt64 = ceil((p_alt - min(p_alt)) / max(p_alt - min(p_alt)) * (cml-1)) + 1;

% Plot things
m_proj('Mercator','long',lonbdy,'lat',latbdy);
m_coast('color','k');
m_states('k');
m_grid('linestyle','none')
title(sprintf('Flight path on %s ', Merge.metadata.date),'fontsize',18)

for a = 1:numel(lat)
    m_line(lon(a), lat(a), 'color',cm((cml+1)-p_alt64(a),:)); % We do the subtraction in the cm index so that high pressures are lower on the colorbar
end

% Now, unfortunately since we've done a custom implementation of coloring,
% we'll need to create our own tick marks. We'll want each tick mark to
% represent the highest pressure that would have that color.

cb = colorbar;
ticks = get(cb,'YTick'); nTicks = numel(ticks);

% We're just going to relabel our existing ticks, so first we find our
% colormap indices that will correspond to those ticks
tick_indices = floor(linspace(1,cml,nTicks));

% Then we want to find the highest pressure that corresponds to each of the
% indices
pressures = cell(size(tick_indices));
for b = 1:numel(tick_indices)
    tmp = p_alt(p_alt64 == tick_indices(b));
    pressures{b} = num2str(round(max(tmp)));
end
pressures = fliplr(pressures); % This will put the highest pressure as the lowest tick on the colorbar
set(cb,'YTickLabel',pressures);
ylabel(cb,'Pressure Altitude (hPa)','fontsize',16)
end

