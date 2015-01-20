%NCAR_LIF_comparison - Plots LIF and NCAR NO2 measurements against each
%other for comparison to examine potential bias.
%
%   The points will be colored by date to check if there is bias on certain
%   days only.

startdate = '9/01/2013';
enddate = '9/30/2013';

LIF_fieldname = 'NO2_MixingRatio_LIF';
NCAR_fieldname = 'NO2_MixingRatio';

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';

DEBUG_LEVEL = 1;

dates = datenum(startdate):datenum(enddate);

no2lif = cell(1,30);
no2ncar = cell(1,30);
datenums = zeros(1,30);
datestrs = cell(1,30);

S=0;
for d=1:numel(dates)
    % Load the merge files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    S=S+1;
    no2lif{S} = remove_merge_fills(Merge,LIF_fieldname);
    no2ncar{S} = remove_merge_fills(Merge,NCAR_fieldname);
    datenums(S) = dates(d); datestrs{S} = datestr(dates(d),2);
end

% Plot coloring by date
figure; hold on
for a=1:S
    datecolor = datenums(a)*ones(size(no2lif{a}));
    scatter(no2lif{a},no2ncar{a},16,datecolor);
end
xlabel('NO2 LIF Mixing Ratio (pptv)','fontsize',14);
ylabel('NO2 NCAR Mixing Ratio (pptv)','fontsize',14);
title(sprintf('LIF vs. NCAR for %s to %s',startdate,enddate),'fontsize',16,'fontweight','bold');
cb = colorbar;
set(cb,'Ticks',datenums(1:S));
set(cb,'TickLabels',datestrs(1:S));

% Concatenate to plot fit line
no2lif_all = cat(2,no2lif{1:S});
no2ncar_all = cat(2,no2ncar{1:S});
nans = isnan(no2lif_all) | isnan(no2ncar_all);
no2lif_all(nans) = []; no2ncar_all(nans) = [];
%[P,R] = polyfit_R2(no2lif_all,no2ncar_all,1);
[slope, int, r] = lsqfitgm(no2lif_all, no2ncar_all);
P(1) = slope; P(2) = int; R = r^2;
xt = get(gca,'xtick');
lfit = line(xt,polyval(P,xt),'linewidth',2,'linestyle','--','color','k');
l1to1 = line(xt,xt,'linewidth',2,'linestyle',':','color','r');
legtext = cell(1,2);
legtext{1} = sprintf('Fit %.4fx + %.2f \n R^2 = %.4f',P(1),P(2),R);
legtext{2} = '1:1';
legend([lfit;l1to1],legtext,'fontsize',12);