%% Create plots for ICE-AGE app
% Run create_satimage_gifs.m to create gif timeseries of Sentinel-2 images
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/melange-characterization-code/', ...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/');

% Plot example data for ASG, ZIM, HLG, and SEK
sites = [{'ILS'},{'ZIM'},{'HLG'},{'SEK'}];
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/';
size_dir = [root_dir,'melange/'];
melt_dir = [root_dir,'iceberg-melt-data/site-csvs/'];
fig_dir = [root_dir,'figures/ICE-AGE/'];

% specify variables (from extract_automated_iceberg_DEM_distributions.m)
month_cmap = cmocean('phase',12);

% Iterate through the sites, loading and plotting the data
for j = 1:length(sites)
    iceberg_fig = figure; set(iceberg_fig,'position',[1050 50 1200 600]);

    % Open the iceberg size distribution timeseries file
    cd([size_dir,char(sites(j))]);
    bergs = readtable([char(sites(j)),'-iceberg-distribution-timeseries.csv']);

    % Plot the size distributions
    subplot(1,2,1);
    Asurf = table2array(bergs(:,1)); Asurf_bin = table2array(bergs(:,2)); 
    for k = 1:size(bergs,2)-2
        %extract the observation date
        header_title = split(char(bergs.Properties.VariableNames(k+2)),'_');
        obs_date = char(header_title(2));
        yr(k) = str2num(obs_date(1:4)); mo(k) = str2num(obs_date(5:6)); 

        %normalize by total melange area
        iceberg_area = Asurf.*table2array(bergs(:,k+2));
        berg_norm = table2array(bergs(:,k+2))./sum(iceberg_area); berg_norm(berg_norm==0) = NaN;

        %plot the data
        if k == 1
            for l = 1:12
                leg_pl(l) = loglog(Asurf,berg_norm,'.-','linewidth',2,'color',month_cmap(l,:)); hold on;
            end
        end
        loglog(Asurf,berg_norm,'.-','linewidth',2,'color',month_cmap(mo(k),:)); hold on;

        clear header_title obs_date iceberg_area berg_norm;
    end
    clear yr mo;
    set(gca,'ylim',[10^-12, 0.05],'xlim',[40 10^6.5],'fontsize',20); grid on;
    xlabel('Iceberg area (m^2)','fontsize',20);
    ylabel('Normalized iceberg count','fontsize',20);
    leg = legend(leg_pl,'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
    set(leg,'location','northeast','fontsize',16);
    
    % Load the melt rate data & plot it
    cd(melt_dir);
    melt = readtable([char(sites(j)),'_basic_iceberg_meltinfo.csv']);
    meltrates = melt.MeltRate_mD__1_; drafts = melt.DraftMean_m_;
    subplot(1,2,2);
    %Find the middate between the start and end dates and use that as the color
    for k = 1:size(melt,1)
        %calculate the mid-date
        start_date = table2array(melt(k,2)); end_date = table2array(melt(k,3)); 
        start_datetime = datetime(num2str(start_date),'InputFormat','yyyyMMdd');
        end_datetime = datetime(num2str(end_date),'InputFormat','yyyyMMdd');
        mid_datetime = mean([start_datetime, end_datetime]);

        %plot the data
        plot(drafts(k),meltrates(k),'sk',...
            'MarkerFaceColor',month_cmap(month(mid_datetime),:),'MarkerEdgeColor',month_cmap(month(mid_datetime),:)); hold on;

        clear start_date* end_date* mid_datetime;
    end
    clear meltrates drafts;
    set(gca,'ylim',[0 1.25],'xlim',[0 400],'fontsize',20); grid on;
    xlabel('Iceberg draft (m)','fontsize',20);
    ylabel('Melt rate (m/d)','fontsize',20);

    % Save the figure
    exportgraphics(iceberg_fig,[fig_dir,char(sites(j)),'-iceberg-sizes-meltrates_scatterplots.tif'],Resolution=600);
    close all;
end
