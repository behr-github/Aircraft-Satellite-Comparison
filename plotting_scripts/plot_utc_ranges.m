%plot_utc_ranges
%
%   Plots ranges specified in a mat file containing a Ranges structure
%   against a field in a Merge data structure by day.  Will also
%   (optionally) plot a map of the flight path with the ranges highlighted
%   red.
%
%   Josh Laughner <joshlaugh5@gmail.com> 24 Jul 2014

range_file = '/Volumes/share/GROUP/INTEX-B/INTEXB_Profile_UTC_Ranges_Inclusive.mat';
merge_dir = '/Volumes/share/GROUP/INTEX-B/Matlab files/';

plot_field = 'ALTITUDE_GPS';

plot_map = true; % set to true to plot a map of the flight path as well
latlim = [10 30]; % Leave these empty to automatically choose lat/lon edges for the map
lonlim = [-120, -60]; 

load coast;
coastlat = lat; coastlon = long; clear lat long

load(range_file);
load('blue_red_cmap');



for a=1:numel(Ranges);
    curr_date = datestr(Ranges(a).Date,29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    fstring = sprintf('*%s_%s_%s.mat',year,month,day);
    fpath = fullfile(merge_dir,fstring);
    files = dir(fpath);
    load(fullfile(merge_dir,files(1).name));%,'Merge'));
    
    data = eval(sprintf('Merge.Data.%s.Values;',plot_field));
    utc = Merge.Data.UTC.Values;
    data_fills = eval(sprintf('Merge.Data.%s.Fill;',plot_field));
    data(data==data_fills) = NaN;
    range_fig = figure; 
    plot(utc,data,'color','k');
    hold on
    yl = get(gca,'ylim');
    colortrack = false(size(utc));
    hold on
    for b=1:size(Ranges(a).Ranges,1)
        r = Ranges(a).Ranges(b,:);
        xx = utc > r(1) & utc < r(2); colortrack(xx) = 1;
        fill([r(1), r(1), r(2), r(2)],[yl(1), yl(2), yl(2), yl(1)],'r','FaceAlpha',0.4)
    end
    title(curr_date,'fontsize',18);
    hold on
    
    if plot_map
        lon = Merge.Data.LONGITUDE.Values - 360;
        fills = Merge.Data.LONGITUDE.Fill;
        lat = Merge.Data.LATITUDE.Values;
        lon(lon==fills)=NaN; lat(lat==fills)=NaN;
        map_fig = figure;
        scatter(lon,lat,16,colortrack); colormap(blue_red_cmap);
        line(coastlon, coastlat,'color','k');
        if isempty(latlim)
            xlim([floor(min(lon(:))/30)*30, ceil(max(lon(:))/30)*30]);
            ylim([floor(min(lat(:))/30)*30, ceil(max(lat(:))/30)*30]);
        else
            xlim(lonlim); ylim(latlim);
        end
        
        mappos = get(map_fig,'position');
        rangepos = get(range_fig,'position');
        
        set(range_fig,'position',[0, rangepos(2:4)]);
        set(map_fig,'position',[rangepos(3)+64,mappos(2:4)]);
        
        title(curr_date,'fontsize',18);
    end
    
    pause;
    close(range_fig);
    if plot_map; close(map_fig); end
    
end