function [ ComparisonSorted ] = make_resid_vs_aod_plots( Comparison, campaign_name )
%make_resid_vs_aod_plots Function to ease the creation of residual vs AOD plots
%   I've been making lots of plots of the difference in satellite and
%   aircraft VCDs vs. AOD.  This function will take a Comparison structure
%   made by Run_all_aer_categories and parse it into a form where all the
%   relevant data is easier to get a hold of.
%
%   If you include a second argument that is a string describing what
%   campaign this is for (does not need to be recognizable by
%   merge_field_names, it's only to go in plot titles) this will also go
%   ahead and make those plots.

ComparisonSorted.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
ComparisonSorted.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

ComparisonSorted.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
ComparisonSorted.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
ComparisonSorted.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
ComparisonSorted.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
ComparisonSorted.coinc.reldiff = (ComparisonSorted.coinc.behrno2 - ComparisonSorted.coinc.airno2)./ComparisonSorted.coinc.airno2 * 100;

ComparisonSorted.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
ComparisonSorted.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
ComparisonSorted.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
ComparisonSorted.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
ComparisonSorted.aerabove.reldiff = (ComparisonSorted.aerabove.behrno2 - ComparisonSorted.aerabove.airno2)./ComparisonSorted.aerabove.airno2 * 100;

ComparisonSorted.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
ComparisonSorted.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
ComparisonSorted.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
ComparisonSorted.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
ComparisonSorted.no2above.reldiff = (ComparisonSorted.no2above.behrno2 - ComparisonSorted.no2above.airno2)./ComparisonSorted.no2above.airno2 * 100;

if nargin > 1
    make_plots(ComparisonSorted, campaign_name);
end

end

function make_plots(ComparisonSorted, campaign_name)
airval_c = ComparisonSorted.coinc.reldiff;
airval_a = ComparisonSorted.aerabove.reldiff;
airval_n = ComparisonSorted.no2above.reldiff;

% Look for fill values and remove
xx = ComparisonSorted.coinc.aod < -1;

if sum(~xx) > 0
    figure;
    xl = [min(min(ComparisonSorted.coinc.aod(~xx)),0), max(ComparisonSorted.coinc.aod(~xx))];
    if isnan(xl(2)); xl(2) = 1; end

    scatter(ComparisonSorted.coinc.aod(~xx),airval_c(~xx));
    set(gca,'xlim',xl)
    oylim = max(abs(get(gca,'ylim')));
    set(gca,'ylim',[-oylim,oylim]);
    set(gca,'fontsize',16)
    [xline,yline,str] = calc_fit_line(ComparisonSorted.coinc.aod(~xx),airval_c(~xx),'regression','rma');
    l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
    legend(l,{str});
    xlabel('AOD');
    ylabel('Percent difference in column');
    title(sprintf('%s: Coincident layers',upper(campaign_name)));
end


% Look for fill values and remove
xx = ComparisonSorted.aerabove.aod < -1;

if sum(~xx) > 0
    xl = [min(min(ComparisonSorted.aerabove.aod(~xx)),0), max(ComparisonSorted.aerabove.aod(~xx))];
    if isnan(xl(2)); xl(2) = 1; end

    scatter(ComparisonSorted.aerabove.aod(~xx),airval_a(~xx));
    set(gca,'xlim',xl);
    oylim = max(abs(get(gca,'ylim')));
    set(gca,'ylim',[-oylim,oylim]);
    set(gca,'fontsize',16)
    [xline,yline,str] = calc_fit_line(ComparisonSorted.aerabove.aod(~xx),airval_a(~xx),'regression','rma');
    l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
    legend(l,{str});
    xlabel('AOD');
    ylabel('Percent difference in column');
    title(sprintf('%s: Aerosol layer above',upper(campaign_name)));
end

% Look for fill values and remove
xx = ComparisonSorted.no2above.aod < -1;

if sum(~xx)
    figure;
    xl = [min(min(ComparisonSorted.no2above.aod(~xx)),0), max(ComparisonSorted.no2above.aod(~xx))];
    if isnan(xl(2)); xl(2) = 1; end

    scatter(ComparisonSorted.no2above.aod(~xx),airval_n(~xx));
    set(gca,'xlim',xl);
    oylim = max(abs(get(gca,'ylim')));
    set(gca,'ylim',[-oylim,oylim]);
    set(gca,'fontsize',16)
    [xline,yline,str] = calc_fit_line(ComparisonSorted.no2above.aod(~xx),airval_n(~xx),'regression','rma');
    l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
    legend(l,{str});
    xlabel('AOD');
    ylabel('Percent difference in column');
    title(sprintf('%s: NO_2 layer above',upper(campaign_name)));
end
end
