function [  ] = examine_aerosol_categorization( profstr )
%examine_aerosol_categorization Plots the LIF and NCAR NO2 profiles against
%aerosol extinction profiles with the 90% column heights marked.

mergeno2dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
mergeaerdir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aerosol_Merges/Properties/';
mergepat = 'Baltimore*%s_%s_%s.mat';
profilefield = 'ProfileSequenceNum';
no2field = 'NO2_LIF';
no2ncarfield = 'NO2_NCAR';
aerosol_field = 'EXTamb532';
aeralt_field = 'GPS_ALT';

dates = datenum('7/1/2011'):datenum('7/31/2011');
for d=1:numel(dates)
    SNO2 = wildcard_load(mergeno2dir,mergepat,dates(d),'Merge');
    SAer = wildcard_load(mergeaerdir,mergepat,dates(d),'Merge');
    if isempty(SNO2) || isempty(SAer); continue; end
    MergeNO2 = SNO2.Merge;
    MergeAer = SAer.Merge;
    
    [no2, utc, presNO2] = remove_merge_fills(MergeNO2,no2field);
    no2ncar = remove_merge_fills(MergeNO2, no2ncarfield);
    [ext,~,presAer] = remove_merge_fills(MergeAer,aerosol_field,'alt',aeralt_field);
    presAer = presAer * 12 * 2.54 * 1e-5;
    xx = utc > local2utc('12:00','est') & utc < local2utc('15:00','est');
    profnums = unique(MergeNO2.Data.(profilefield).Values(xx));
    profnums(profnums<1) = [];
    
    %fieldid = 'CoincidentHigh';
    %alts = profstr.([fieldid,'Alt']);
    for pn = 1:numel(profnums);
        p = profnums(pn);
        %if all(p~=profstr.(fieldid)); continue; end 
        %ff = p == alts(:,1);
        ppa = MergeAer.Data.(profilefield).Values == p;
        ppn = MergeNO2.Data.(profilefield).Values == p;
        [no2bins, no2alts] = bin_rolling_vertical_profile(presNO2(ppn),no2(ppn),0.5,0.1);
        [ncarbins, ncaralts] = bin_rolling_vertical_profile(presNO2(ppn), no2ncar(ppn), 0.5, 0.1);
        [aerbins, airalts] = bin_rolling_vertical_profile(presAer(ppa),ext(ppa),0.5,0.1);
        figure;
        hax = plotxx(aerbins, airalts, no2bins, no2alts, {'Aerosol extinction at 532 nm','NO2 pptv'},{'Alt (km)','Alt (km)'});

        oldlim = get(hax(1),'xlim'); 
        newlim = [min(0,min(aerbins)),max(0,max(oldlim))];
        set(hax(1),'xlim', newlim);
        %line([0, oldlim(2)/2, oldlim(2)],[alts(ff,3),alts(ff,3),alts(ff,3)],'linestyle','--','color','k','linewidth',2,'marker','x','parent',hax(1)); 
        
        oldlim = get(hax(2),'xlim'); 
        newlim = [min(0,min(no2bins)),max(0,max(oldlim))];
        set(hax(2), 'xlim', newlim);
        line(ncarbins, ncaralts, 'color', [0 0.7 0], 'parent',hax(2));
        %line([0, oldlim(2)],[alts(ff,2),alts(ff,2)],'linestyle','--','color','r','linewidth',2,'parent',hax(2));
        
        ylims1 = get(hax(1),'ylim');
        ylims2 = get(hax(2),'ylim');
        newylim = [0, max([ylims1, ylims2])];
        set(hax(1),'ylim',newylim);
        set(hax(2),'ylim',newylim);
        
        dt = MergeNO2.metadata.date;
        title(sprintf('%s Prof Num %d',dt,p));
    end
    if ~isempty(get(0,'children'))
        tilefigs;
        pause;
        close all
    end
end


end

