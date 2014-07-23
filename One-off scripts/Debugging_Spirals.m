%Debugging_Spirals
%
%   Sequentially creates maps with BEHR NO2 maps, overlaid with the flight
%   path for that day and the pixel boundaries, plus a second figure where
%   the pixel centers are marked with the difference in columns

date_start = '07/01/2011';
date_end = '07/31/2011';
no2field = 'NO2_LIF';   %For Baltimore, this is NO2_LIF or NO2_NCAR
%For CA/TX, this is NO2_MixingRatio_LIF or NO2_MixingRatio

tz = 'est';

states={'md'};

% Which plots to show
avgNO2 = true;
deltaNO2 = true;
dailySatNO2 = true;
stratNO2 = false;

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR/';

behr_map_file = '/Users/Josh/Documents/MATLAB/Figures/Sat Verification/Debugging Spiral Method/BEHR Columns July 2011.mat';

DEBUG_LEVEL = 1;

dates = datenum(date_start):datenum(date_end);

load('blue_red_cmap.mat');
load(behr_map_file);
latbdy = [min(LAT_GRID(:)), max(LAT_GRID(:))];
lonbdy = [min(LON_GRID(:)), max(LON_GRID(:))];
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('OMI_BEHR_*%s%s%s.mat',year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    behr_files = dir(fullfile(behr_dir,behr_filename));
    if numel(behr_files)==1
        load(fullfile(behr_dir,behr_files(1).name),'Data')
    elseif isempty(behr_files)
        if DEBUG_LEVEL > 1; fprintf('No BEHR file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of BEHR files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    S=0;
    lon = cell(1,4); lat = cell(1,4); omino2 = cell(1,4); behrno2 = cell(1,4); airno2 = cell(1,4);
    airno2err = cell(1,4); cov = cell(1,4); quality = cell(1,4); db = struct;
    for swath=1:numel(Data)
        S=S+1;
        [lon{S}, lat{S}, omino2{S}, behrno2{S}, airno2{S}, airno2err{S}, cov{S}, quality{S}, db(S).db] = spiral_verification(Merge,Data(swath),tz,'DEBUG_LEVEL',DEBUG_LEVEL,'no2field',no2field);
    end
    
    lonall = cat(1,lon{:}); latall = cat(1,lat{:});
    if isempty(lonall)
        continue
    else
        
        % The satellite map
        if avgNO2
            fmap = figure; m_proj('Mercator','long',lonbdy,'lat',latbdy);
            m_pcolor(LON_GRID, LAT_GRID, NO2_GRID); shading flat
            m_states('w'); colorbar; caxis([0, 1e16]);
            m_grid('linestyle','none');
            
            for a=1:numel(db)
                s = size(db(a).db.loncorn);
                if s(2) > 0;
                    for b=1:s(2)
                        m_line([db(a).db.loncorn(:,b); db(a).db.loncorn(1,b)], [db(a).db.latcorn(:,b); db(a).db.latcorn(1,b)],'color','k','linewidth',2)
                    end
                end
            end
            
            m_line(Merge.Data.LONGITUDE.Values-360, Merge.Data.LATITUDE.Values, 'color', 'r', 'linewidth',1);
            title(sprintf('%s',Merge.metadata.date),'fontsize',20)
            oldposmap = get(fmap,'Position'); set(fmap,'Position',[1,1, oldposmap(3:4)*1.5]);
        end
        
        % The difference map
        if deltaNO2
            fcent = state_outlines(states{:});
            delta = cat(1,behrno2{:}) - cat(1,airno2{:});
            for a=1:numel(db)
                s = size(db(a).db.loncorn);
                if s(2) > 0;
                    for b=1:s(2)
                        line([db(a).db.loncorn(:,b); db(a).db.loncorn(1,b)], [db(a).db.latcorn(:,b); db(a).db.latcorn(1,b)],'color','b','linewidth',2)
                    end
                end
            end
            hold on;
            scatter(lonall, latall, 20, delta);
            dx = abs(lonbdy(2)-lonbdy(1))/8;
            for a=1:numel(delta)
                % Find if any previous texts are at the same place
                if a>1
                    samept = sum(abs(lonall(1:a-1)-lonall(a)) < 0.1 & abs(latall(1:a-1)-latall(a)) < 0.1);
                else
                    samept = 0;
                end
                text(lonall(a)+dx/2+samept*dx, latall(a), sprintf('%+.1e',delta(a)),'BackgroundColor',[0.6 0.6 0.6]);
            end
            xlim(lonbdy); ylim(latbdy);
            colorbar; colormap(blue_red_cmap); caxis([-5e15 5e15]);
            title(sprintf('%s',Merge.metadata.date),'fontsize',20)
            oldposcent = get(fcent,'Position');
            set(fcent,'Position',[oldposmap(3)*2,1,oldposcent(3:4)*1.5]);
        end
        
        %A satellite map for that day
        if dailySatNO2
            cb = no2_column_map_2014(curr_date, curr_date, lonbdy, latbdy, 'rowanomaly', 'AlwaysByRow', 'projection', 'Mercator',...
                'behrdir',behr_dir,'fileprefix','OMI_BEHR_omiCloudAMF_','cbrange',[0 1e16]);
            fsatday = gcf;
        end
        
        %Plot the stratospheric columns
        if stratNO2
            fstrat = state_outlines(states{:});
            SNO2 = [];
            for a=1:numel(db)
                s = size(db(a).db.loncorn);
                if s(2) > 0;
                    SNO2 = cat(1,SNO2,db(a).db.strat_NO2);
                    for b=1:s(2)
                        line([db(a).db.loncorn(:,b); db(a).db.loncorn(1,b)], [db(a).db.latcorn(:,b); db(a).db.latcorn(1,b)],'color','b','linewidth',2)
                    end
                end
            end
            hold on;
            scatter(cat(1,lon{:}), cat(1,lat{:}), 20, SNO2); colorbar;
            xlim(lonbdy); ylim(latbdy);
            title(sprintf('%s Stratospheric NO2',Merge.metadata.date),'fontsize',20)
        end
        
        pause;
        if deltaNO2; close(fcent); end
        if avgNO2; close(fmap); end
        if dailySatNO2; close(fsatday); end
        if stratNO2; close(fstrat); end
    end
    
end