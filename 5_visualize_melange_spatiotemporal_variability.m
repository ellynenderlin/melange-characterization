%%% Create plots & maps of GrIS melange data

clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/miscellaneous/melange-characterization-code/');

%specify root path
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
% root_dir = '/Volumes/Jokulhaup_5T/Greenland-melange/';

%define custom parameters for size distributions
transect_inc = 1000; %distance between transects along the centerline (meters)
zthresh = 3; %cutoff elevation (m)

%Thickness parameters:
zcutoff = zthresh; %elevation threshold below which to ignore icebergs (m)
rho_i = 900; %ice density (kg/m^3)
rho_sw = 1026; %water density (kg/m^3)
Hcutoff = round((rho_sw/(rho_sw-rho_i))*zcutoff); %H threshold for figure naming

%Size distribution parameters
nthresh = 1e-6; % set small number bin cutoff (n1 must be greater than this value)
ARcomp.best.autoALL = 2; % iceberg aspect ratio 
vthresh = (1/4)*pi*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*zthresh).^2; %don't include this bin size or smaller in curve fitting
dplawthresh = 10^5; % upper bound on the intercept for the dummy powerlaw
norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'
normalize_exp = 1.5; % Increase to weight residuals towards end of the curve, minimum = 1

%Velocity parameters: ideally only use velocities from time periods with
%DEMs to create the seasonal profiles but you can also use all velocities
%from that season (even if there is no DEM for that season-year combination)
vfilter = 'annual'; %options: 'annual' OR 'all'
vdtmin = 15; %maximum image acquisition separation for velocities
vdtmax = 60; %maximum image acquisition separation for velocities
dvterm = 'yes'; %include velocities from the terminus

%decide whether to sample data at (1) fixed locations or (2) relative
%locations with respect to the terminus positions in the DEMs
sampling = 'dated'; % options: 'dated' OR 'fixed';

%customize visualization
years = 2011:1:2023; yr_cmap = cmocean('matter',length(years)+1); yr_cmap = yr_cmap(2:end,:);
seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MAM','JJA','SON'};
mo_cmap = cmocean('phase',12); 
% seas_cmap = [194,165,207; 123,50,148; 0,136,55; 166,219,160]/255; %purple/green
seas_cmap = [5,113,176; 202,0,32; 244,165,130; 146,197,222]/255; %blue/orange/red
%set up site-aggregated plots to be arranged somewhat geographically but with a hole for the site map to be inserted
% rows = 9; cols = 2; %9 sites in west Greenland
% plot_locs = [1,3,5,7,9,11,13,15,17,18,16,14,12,2];
rows = 7; cols = 3; %9 sites in west Greenland
plot_locs = [2,1,4,7,10,13,16,19,20,21,18,15,12,3];
close all;

%identify the site folders
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 3
        sitenames = [sitenames; sites(i).name];
    end
end
%if you have a preferred order for the sites in plots, specify it here
geo_order = [{'ULS'},{'KOG'},{'ASG'},{'ILS'},{'UNS'},{'SAS'},{'UMS'},...
    {'KGS'},{'SEK'},{'HLG'},{'MGG'},{'MDG'},{'DJG'},{'ZIM'}];
geo_names = [{'Ullip'},{'Nuussuup'},{'Nunatakassaap'},...
    {'Illullip'},{'Upernavik North'},{'Salliarutsip'},{'Umiammakku'},...
    {'Kangilliup'},{'Sermeq Kujalleq'},{'Helheim'},{'Nigertiip Apusiia'},...
    {'Magga Dan'},{'Daugaard Jensen'},{'Zachariae Isstrom'}];
for j = 1:length(geo_order)
    geo_ind(j) = find(contains(string(sitenames),geo_order(j)));
end
big3 = [{'SEK'},{'HLG'},{'ZIM'}];

alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m',...
    'n','o','p','q','r','s','t','u','v','w','x','y','z'];


%% load the saved data compiled for all study sites
load([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat']);
for j = 1:length(MP)
    saved_sites(j,:) = MP(j).name;
end

%% create overview plots of melange characteristics
close all; drawnow;

% %load the melange characteristic dataset created by compile_melange_characteristics.m
% load([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat']);
% for j = 1:length(MP)
%     saved_sites(j,:) = MP(j).name;
% end

% create summary figures
summaryfig = figure; set(summaryfig,'position',[-1650 50 900 1200]);
for j = 1:length(geo_ind)
    figure(summaryfig);
    eval(['subH',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
    % figure(Vfig);
    % eval(['subV',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
end
missfig = figure; set(missfig,'position',[1050 50 900 500]); 
% subm1 = subplot(2,1,1); subm2 = subplot(2,1,2);

% %If adding new dates, profiles need to be removed from the structure
% %because the inland reference may have changed
% %(will throw errors when concatenating)
% for j = 1:length(MP)
%     MP(j).Z = rmfield(MP(j).Z,{'dist','Hseas'});
%     MP(j).V = rmfield(MP(j).V,{'dist','Vseas','dVdx'});
% end

%create dummy vectors to hold all buttressing data for histogram plotting
BM_annual = []; BM_character = []; BA_annual = []; BA_character = [];
%iterate
disp(['Creating profile plots ignoring icebergs thinner than ',num2str(Hcutoff),'m & using ',vfilter,' speeds w/ dt ~',num2str(round([vdtmin,vdtmax])),' days']);
for l = 1:length(geo_names)
    %identify the index in the structure for the specified sitename
    j = find(strcmp(string(saved_sites),string(geo_order(l))));

    %specify site plotting conventions
    site_abbrev = MP(j).name; disp(site_abbrev);
    seasons = MP(j).D.months;
    disp(site_abbrev);
    transect_inc = round(nanmean(sqrt((diff(MP(j).V.X).^2 + diff(MP(j).V.Y).^2)))/100)*100;

    %isolate the dates for data sorting
    D_dir = [root_dir,site_abbrev,'/'];
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
    end
    [berg_dates,Havg,packing,inland_idx,seaward_idx,seaward_ext,H_seas,pack_seas] = extract_seasonal_binned_DEM_profiles(D_dir,site_abbrev,term_trace,MP(j).Z.termX,MP(j).Z.termY,MP(j).V.X,MP(j).V.Y,ARcomp.best.autoALL,zcutoff,sampling,years,seasons);
    %create derived variables from elevation-based profile outputs
    seaward_meanidx = round(mean(seaward_ext));
    for p = 1:size(berg_dates,2)
        zyrs(p) = year(berg_dates(p)); zmos(p) = month(berg_dates(p));
    end
    Hdist = MP(j).Z.dist; 

    %if using all velocities for seasonal profiles, not just from DEM
    %years, create a timeseries of the inland_idx based off the terminus
    %delineations
    for k = 1:length(MP(j).T.termX)
        try 
            Tdatetime(k) = datetime(MP(j).T.date{k},'InputFormat','yyyyMMdd');
        catch
            Tdatetime(k) = datetime(MP(j).T.date{k},'InputFormat','yyyyMMdd');
        end
        Tdate(k) = convert_to_decimaldate(char(MP(j).T.date{k}));
        Tdatetime(k) = datetime(MP(j).T.date{k},'InputFormat','yyyyMMdd');
        Tyrs(k) = year(Tdatetime(k)); Tmos(k) = month(Tdatetime(k));
        dists = sqrt((MP(j).T.termX(k)-MP(j).V.X).^2 + (MP(j).T.termY(k)-MP(j).V.Y).^2);
        [~,ia] = sort(dists);
        term_idx(k) = max(ia(1:2));
        clear dists ia;
    end
    for p = 1:length(years)
        yr_idx = find(Tyrs == years(p));
        vterm_idx(p) = round(median(term_idx(yr_idx)));
    end

    %grab velocity timeseries for the transect-centerline intersection
    %points & create seasonally-averaged climatologies
    if contains(dvterm,'y')
        inland_vel = max(inland_idx); %grab data from the terminus
        ref_adjust = 1; %shift the referencing for the relative coordinate system
    else
        inland_vel = max(inland_idx)-1; %only get speeds from the melange
        ref_adjust = 0;
    end
    vel_dir = [root_dir,site_abbrev,'/velocities/'];
    vel_seas = create_seasonal_speed_profiles(vel_dir,inland_idx,inland_vel,ref_adjust,vterm_idx,term_idx,Tyrs,Tmos,vdtmin,vdtmax,vfilter,years,seasons,zyrs,zmos);
    vdist = MP(j).V.dist;

    %loop through the subsetted size distributions and compute seasonal
    %averages nearest the terminus & a fixed relative distance down-fjord
    Dsubs = dir([root_dir,site_abbrev,'/',site_abbrev,'*-iceberg-distribution-subsets.csv']);
    for p = 1:length(Dsubs)
        D = readtable([root_dir,site_abbrev,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        berg_datestrings(p,:) = Dsubs(p).name(length(site_abbrev)+2:length(site_abbrev)+9);
        byrs(p) = year(berg_dates(p)); bmos(p) = month(berg_dates(p));
        if p == 1 %grab the bin size info from the first file (same in all files)
            bHs = ((2/ARcomp.best.autoALL)*sqrt(table2array(D(:,1))/pi())); 
            bhs = ((rho_sw-rho_i)/rho_sw)*bHs; dsize = round(nanmean(diff(bhs)));
            berg_maxz = bhs+0.5*dsize; size_ind = find(berg_maxz <= zcutoff,1,'last');
        end

        if ~isnan(inland_idx(p))
            %full size distributions at two points
            bergdist_inland(p,:) = table2array(D(:,(inland_idx(p)-1)+2))'; %near-terminus
            bergdist_seaward_setdx(p,:) = table2array(D(:,seaward_meanidx+2))'; %standard distance from near-terminus
        else
            bergdist_inland(p,:) = NaN(size(table2array(D(:,1))))'; %near-terminus
            bergdist_seaward_setdx(p,:) = NaN(size(table2array(D(:,1))))'; %standard distance from near-terminus
        end
        bergdist_seaward_end(p,:) = table2array(D(:,seaward_ext(p)+2))'; %seaward melange limit in DEM
        clear D;
    end
    bergdist_inland(bergdist_inland==0) = NaN; bergdist_inland_norm = bergdist_inland./sum((bergdist_inland.*MP(j).D.area'),2,'omitnan');
    bergdist_seaward_setdx(bergdist_seaward_setdx==0) = NaN; bergdist_seaward_setdx_norm = bergdist_seaward_setdx./sum((bergdist_seaward_setdx.*MP(j).D.area'),2,'omitnan');
    bergdist_seaward_end(bergdist_seaward_end==0) = NaN; bergdist_seaward_end_norm = bergdist_seaward_end./sum((bergdist_seaward_end.*MP(j).D.area'),2,'omitnan');
    %create seasonal average distributions for each fjord bin
    % for k = 1:4
    %     bergdist_seas(k,:,1) = nanmean(bergdist_seaward_setdx_norm(ismember(bmos,seasons(k,:))==1,:),1); %could change to bergdist_seaward_end_norm for moving end
    %     bergdist_seas(k,:,2) = nanmean(bergdist_inland_norm(ismember(bmos,seasons(k,:))==1,:),1);
    % end
    %dummy matrices to hold seasonal average distributions for each year
    bergdist_inland_seas = NaN(4,size(bergdist_inland,2),length(years));
    bergdist_seaward_seas = NaN(4,size(bergdist_inland,2),length(years));
    for p = 1:length(years)
        yr_idx = find(byrs == years(p));
        if ~isempty(yr_idx)

            %isolate the distributions for that year
            for k = 1:length(yr_idx)
                bi_temp(k,:) = bergdist_inland_norm(yr_idx(k),:);
                bs_temp(k,:) = bergdist_seaward_setdx_norm(yr_idx(k),:); %could change to bergdist_seaward_end_norm for moving end
            end

            %calculate seasonal average
            mos_yr = bmos(yr_idx);
            for k = 1:4
                if ~isempty(ismember(mos_yr,seasons(k,:)))
                    bergdist_inland_seas(k,:,p) = nanmean(bi_temp(ismember(mos_yr,seasons(k,:))==1,:),1); %inland
                    bergdist_seaward_seas(k,:,p) = nanmean(bs_temp(ismember(mos_yr,seasons(k,:))==1,:),1); %seaward
                end
            end

            clear mos_yr b*_temp;
        end
        clear yr_idx
    end

    % %plot the thickness distributions for near the terminus
    % if zcutoff == 3
    %     Hdistfig = figure; set(gcf,'position',[600 600 400 400]);
    %     for p = 1:size(bergdist_inland,1)
    %     plot(bHs(size_ind:end),bergdist_inland(p,size_ind:end),'-k','linewidth',1.5); hold on;
    %     end
    % end
    % figure(Hdistfig);
    % disp('Check the thickness distributions make sense, then close the figure to advance'); disp(' ');
    % uiwait %advance only after figure is closed

    %average the seasonal average datasets
    bergdist_seas = NaN(4,size(bergdist_inland,2),2);
    for k = 1:4
        if ~isempty(ismember(bmos,seasons(k,:)))
            bergdist_seas(k,:,1) = nanmean(bergdist_inland_seas(k,:,:),3); %inland
            bergdist_seas(k,:,2) = nanmean(bergdist_seaward_seas(k,:,:),3); %seaward
        end
    end
    disp('extracted distributions');

    %load iceberg distribution profiles
    D = readtable([root_dir,site_abbrev,'/models/',site_abbrev,'-powerlaw-slope_',num2str(zcutoff),'m-zcutoff_seasonal-profiles.csv'],"VariableNamingRule","preserve");
    size_plslope = table2array(D(:,2:end)); clear D;
    D = readtable([root_dir,site_abbrev,'/models/',site_abbrev,'-powerlaw-bergybit-mispredict_',num2str(zcutoff),'m-zcutoff_seasonal-profiles.csv'],"VariableNamingRule","preserve");
    res_A = table2array(D(:,2:end)); clear D;
    %find the last binned seasonal distribution with data
    for k = 1:4
        if sum(~isnan(res_A(k,:))) > 0
            res_Alim(k) = find(~isnan(res_A(k,:))==1,1,'last');
        else
            res_Alim(k) = NaN;
        end
    end
    
    %Site-specific plots
    sitefig = figure; set(sitefig,'position',[50 850 1200 1000]);
    % subz = subplot(4,2,[1:2]); subv = subplot(4,2,[3:4]); %references for thickness & speed profiles
    subis = subplot(4,2,1); subss = subplot(4,2,2); %inland and seaward size distributions
    subdv = subplot(4,2,[3:4]); subp = subplot(4,2,[5:6]); subm = subplot(4,2,[7:8]); %references for speed & packing density profiles & missfit scatterplot
    max_xlim = max([Hdist(find(sum(sum(~isnan(H_seas),2),3)>0,1,'last'))]);
    plot_ind = find(geo_ind == geo_ind(l));
    % disp([MP(j).name,': j = ',num2str(j),', plot_ind = ',num2str(plot_ind)])

    %add thickness profiles to figures
    for k = 1:4
        if k <= 2
            fill_alpha = 0.2;
        else
            fill_alpha = 0.4;
        end

        %relative to moving terminus
        Hmean = nanmean(H_seas(:,k,:),3);
        % Hmax = (nanmean(H_seas(:,k,:),3)+std(H_seas(:,k,:),0,3,'omitnan')); 
        % Hmin = (nanmean(H_seas(:,k,:),3)-std(H_seas(:,k,:),0,3,'omitnan')); 
        Hmax = max(H_seas(:,k,:),[],3); 
        Hmin = min(H_seas(:,k,:),[],3); 

        %plot
        Hmax_idx = find(~isnan(Hmax)==1); Hmin_idx = find(~isnan(Hmin)==1);
        if sum(~isnan(Hmean)) > 0
            % %site fig
            % figure(sitefig); subplot(subz);
            % % fill([zdist(Hmax_idx), fliplr(zdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
            % %     seas_cmap(k,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            % % pz(k) = plot(zdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            % fill([Hdist(Hmax_idx), fliplr(Hdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
            %     seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            % pz(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            % %add symbols
            % if ismember(MP(j).name,big3)
            %     plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            % else
            %     plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            % end

            %regional fig
            figure(summaryfig);
            if ~isempty(plot_ind)
                %navigate to the subplot
                eval(['subplot(subH',num2str(geo_ind(l)),');']);
                ax1 = gca;

                %plot the thickness profiles
                yyaxis left
                fill([Hdist(Hmax_idx), fliplr(Hdist(Hmin_idx))]',[Hmax(Hmax_idx); flipud(Hmin(Hmin_idx))],...
                        seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                pH(k) = plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
                if ismember(MP(j).name,big3)
                    % overlay diamond symbols
                    plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'d','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                else
                    % overlay square symbols
                    plot(Hdist(~isnan(Hmean))',Hmean(~isnan(Hmean)),'s','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                end
            end
        else
           % figure(sitefig); subplot(subz);
           % pz(k) = plot(NaN,NaN,'-','color',seas_cmap(k,:),'linewidth',3); hold on; 
        end
        H_ylim(k,:) = [min(Hmean), max(Hmean)];
        clear Hmean Hmax* Hmin*;
    end
    %format the thickness subplot in the site figure
    % figure(sitefig); subplot(subz);
    % seas_leg = legend(pz,season_names);
    % set(gca,'fontsize',16); drawnow; grid on; 
    % set(subz,'xlim',[0,max_xlim],...
    %     'xticklabel',[],'ylim',[floor(min(H_ylim(:,1))/10)*10 ceil(max(H_ylim(:,2))/10)*10]); 
    % ylims = get(gca,'ylim'); 
    % if range(ylims) <= 50
    %     set(gca,'ytick',[floor(min(H_ylim(:,1))/10)*10:10:ceil(max(H_ylim(:,2))/10)*10]);
    % else
    %     set(gca,'ytick',[floor(min(H_ylim(:,1))/10)*10:20:ceil(max(H_ylim(:,2))/10)*10]);
    % end
    % clear ylims;
    % ylabel('Thickness (m)','fontsize',16); %xlabel('Distance from terminus (km)','fontsize',16); 
    % subz_pos = get(subz,'position'); 
    % seas_leg.Location = 'northoutside'; seas_leg.Orientation = 'horizontal'; %move the legend
    % set(subz,'position',subz_pos); %resize the subplot
    %format the site subplot in the thickness figure
    if ~isempty(plot_ind)
        figure(summaryfig); eval(['subplot(subH',num2str(geo_ind(l)),');']);
        yyaxis left;
        set(gca,'xlim',[0, 25000],'xticklabel',[],...
            'ylim',[0 ceil(max(H_ylim(:,2))/10)*10],'fontsize',14); drawnow;
        gca_pos = get(gca,'position'); ylims = get(gca,'ylim'); 
        %slightly adjust ytick plotting depending on zcutoff
        if zcutoff == 3
            %show plot like it is split from speed profile plot
            if max(ylims) <= 50
                set(gca,'ylim',[-10,50],'ytick',[20:10:50]);
                plot([0, 25000],[20,20],'-k','linewidth',1.5);
            else
                set(gca,'ylim',[-50,100],'ytick',[25:25:100]);
                plot([0, 25000],[25,25],'-k','linewidth',1.5);
            end
        else
            if max(ylims) <= 30
                set(gca,'ylim',[-10,30],'ytick',[10:10:30]);
                plot([0, 25000],[10,10],'-k','linewidth',1.5);
            elseif max(ylims) <= 50 && max(ylims) > 30
                set(gca,'ylim',[-30,50],'ytick',[20:10:50]);
                plot([0, 25000],[10,10],'-k','linewidth',1.5);
            elseif max(ylims) <= 100 && max(ylims) > 50
                set(gca,'ylim',[-80,100],'ytick',[25:25:100]);
                plot([0, 25000],[10,10],'-k','linewidth',1.5);
            else
                set(gca,'ylim',[-90,150],'ytick',[50:50:150]);
                plot([0, 25000],[30,30],'-k','linewidth',1.5);
            end
        end
        % yticks = get(gca,'ytick'); set(gca,'ytick',yticks(2:end),'yticklabel',yticks(2:end));
        ax1.YAxis(1).Color = 'k'; 
        set(ax1,'box','on'); ax1.LineWidth = 1.5;
        %add labels
        if plot_locs(plot_ind) >= rows*cols - (cols-1) || (plot_locs(plot_ind) >= 2 && plot_locs(plot_ind) <= 3)
            xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/1000); clear xticks;
            % set(gca,'yticklabel',yticks);
            xlabel('Distance from terminus (km)','fontsize',14); 
            if plot_locs(plot_ind) == rows*cols - (cols-1)
                ax1.YAxis(1).Label.String = 'Thickness (m)';
                ax1.YAxis(1).Label.Position = [-2.7066e3 250 -1];
                ax1.YAxis(1).Label.FontSize = 14;
                % ylabel('Thickness (m)','fontsize',14);
            end
        end
        clear ylims yticks;
        %add a legend
        if j == length(MP)
            H_leg = legend(pH,season_names);
            subH_pos = get(gca,'position');
            H_leg.Orientation = 'horizontal'; %move the legend
            H_leg.Position = [0.38 0.24 0.25 0.035]; 
            set(gca,'position',subH_pos,'fontsize',14); %resize the subplot
            drawnow;
        end

        %slightly shift plot locations to fit inset map in middle
        if plot_ind > 1 && plot_ind <= 8
            set(gca,'position',[gca_pos(1)-0.06 gca_pos(2) gca_pos(3)+0.04 0.068]);
        elseif plot_ind > 9
            set(gca,'position',[gca_pos(1)+0 gca_pos(2) gca_pos(3)+0.04 0.068]);
        else
            set(gca,'position',[gca_pos(1)-0.03 gca_pos(2) gca_pos(3)+0.04 0.068]);
        end
    end
    
    %add speed-based profiles to figures
    for k = 1:4
        if k <= 2
            fill_alpha = 0.2;
        else
            fill_alpha = 0.4;
        end

        %relative to moving terminus
        v_mean = nanmean(vel_seas(:,k,:),3)./365;
        % vmax = (nanmean(vel_seas(:,k,:),3)+std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        % vmin = (nanmean(vel_seas(:,k,:),3)-std(vel_seas(:,k,:),0,3,'omitnan'))./365; 
        vmax = max(vel_seas(:,k,:),[],3)./365; 
        vmin = min(vel_seas(:,k,:),[],3)./365; 
        
        %mask out values near zero because they are so much slower than the
        %melange velocities that they must be bad sea ice tracking
        bottomout = find(v_mean < 0.5,1,'first');
        v_mean(bottomout:end) = NaN; vmax(bottomout:end) = NaN; vmin(bottomout:end) = NaN; 

        %plot
        vmax_idx = find(~isnan(vmax)==1); vmin_idx = find(~isnan(vmin)==1);
        dv_idx = find(~isnan(MP(j).V.dVdx(k,:))==1);
        if sum(~isnan(MP(j).V.dVdx(k,:))) ~= 0
            figure(sitefig); 
            % subplot(subv);
            % fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
            %     seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            % pv(k) = plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            % %add symbols
            % if ismember(MP(j).name,big3)
            %     plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            % else
            %     plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            % end
            subplot(subdv);
            pv(k) = plot(Hdist(dv_idx)',MP(j).V.dVdx(k,dv_idx),'-','color',seas_cmap(k,:),'linewidth',3); hold on;
            %add symbols
            if ismember(MP(j).name,big3)
                plot(Hdist(dv_idx)',MP(j).V.dVdx(k,dv_idx),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            else
                plot(Hdist(dv_idx)',MP(j).V.dVdx(k,dv_idx),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on;
            end

            %regional fig
            figure(summaryfig);
            if ~isempty(plot_ind)
                %navigate to the subplot
                eval(['subplot(subH',num2str(geo_ind(l)),');']);

                %plot the speed profiles
                yyaxis right;
                fill([vdist(vmax_idx), fliplr(vdist(vmin_idx))]',[vmax(vmax_idx); flipud(vmin(vmin_idx))],...
                    seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
                if k >= 2 && sum(~isnan(v_mean))>0
                    pV(k-1) = plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
                else
                    plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
                end
                %add symbols
                if ismember(MP(j).name,big3)
                    %overlay diamond symbols
                    plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'d','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                else
                    %overlay square symbols
                    plot(vdist(~isnan(v_mean))',v_mean(~isnan(v_mean)),'s','color',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'markersize',5); hold on;
                end
            end
        else
            % figure(sitefig); subplot(subv);
            figure(sitefig); subplot(subdv);
            pv(k) = plot(NaN,NaN,'-','color',seas_cmap(k,:),'linewidth',3); hold on;
        end
        clear v_mean vmax* vmin*;
    end
    %format the speed subplot in the site figure
    figure(sitefig); 
    % subplot(subv);
    % set(gca,'fontsize',16); grid on; drawnow;
    % set(subv,'xlim',[0,max_xlim]);
    % xticks = get(subv,'xtick'); set(subv,'xticklabel',xticks/1000); clear xticks; 
    % ylims = get(subv,'ylim'); set(subv,'ylim',[0 max(ylims)]); 
    % yticks = get(gca,'ytick'); 
    % if length(yticks) < 3
    %     if range(ylims) <= 20
    %         set(gca,'ytick',[0:5:max(ylims)]);
    %     else
    %         set(gca,'ytick',[0:10:max(ylims)]);
    %     end
    % end
    % clear ylims yticks;
    % xlabel('Distance from terminus (km)','fontsize',16); ylabel('Speed (m/d)','fontsize',16);
    % pos = get(subv,'position'); %set(subv,'position',[pos(1) pos(2)+0.05 pos(3) pos(4)]);
    subplot(subdv);
    seas_leg = legend(pv,season_names);
    set(gca,'fontsize',16); drawnow; grid on; 
    set(subdv,'xlim',[0,max_xlim]);
    xticks = get(subdv,'xtick'); set(subdv,'xticklabel',xticks/1000); clear xticks; 
    ylims = get(subdv,'ylim'); set(subdv,'ylim',[min(ylims) max(ylims)]); 
    yticks = get(gca,'ytick');
    clear ylims;
    ylabel('Strain rate (d^{-1})','fontsize',16); xlabel('Distance from terminus (km)','fontsize',16); 
    pos = get(subdv,'position'); 
    subdv_pos = get(subdv,'position'); 
    seas_leg.Location = 'northoutside'; seas_leg.Orientation = 'horizontal'; %move the legend
    set(subdv,'position',pos); %resize the subplot

    %format the site subplot in the speed figure
    if ~isempty(plot_ind)
        figure(summaryfig); eval(['subplot(subH',num2str(geo_ind(l)),');']); yyaxis right;
        set(gca,'xlim',[0,25000],'xticklabel',[],'fontsize',14); drawnow;
        % gca_pos = get(gca,'position');
        ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]);
        if max(ylims) <= 35
            set(gca,'ylim',[0,70],'ytick',[0:10:30]);
        % elseif max(ylims) <= 40
        %     set(gca,'ylim',[0,40],'ytick',[0:20:40]);
        else
            set(gca,'ylim',[0,140],'ytick',[0:20:60]);
        end
        yticks = get(gca,'ytick'); 
        ax1.YAxis(2).Color = 'k'; 
        if plot_locs(plot_ind) >= rows*cols - (cols-1)
            xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/1000); clear xticks;
            set(gca,'yticklabel',yticks);
            xlabel('Distance from terminus (km)','fontsize',14);
            if plot_locs(plot_ind) == rows*cols
                ax1.YAxis(2).Label.String = 'Speed (m/d)';
                ax1.YAxis(2).Label.Position = [2.7559e4 550 -1];
                ax1.YAxis(2).Label.FontSize = 14;
                % ylabel('Speed (m/d)','fontsize',14);
            end
        end
        clear ylims yticks;
    end
    gca_pos = get(gca,'position');
    set(gca,'position',[gca_pos(1) gca_pos(2) gca_pos(3) 1.4*gca_pos(4)]);
    drawnow; clear plot_ind pv pV;
    
    %add subplots that generally show the location of the size distributions
    %inland bin
    figure(sitefig); %subplot(4,2,5);
    subplot(subis);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,1),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',16); grid on;
    xlabel('Surface area (m^2)','fontsize',16); ylabel('Iceberg count','fontsize',16);
    % text(10000,10^-2,'near-terminus','fontsize',16);
    %seaward bin
    % subplot(4,2,6);
    subplot(subss);
    for k = 1:4
        loglog(MP(j).D.area,bergdist_seas(k,:,2),'-','color',seas_cmap(k,:),'linewidth',2); hold on;
    end
    set(gca,'ylim',[10^-10 10^-1],'xlim',[10^1 10^7],...
        'ytick',[10^-9,10^-5,10^-1],...
        'xtick',[10^2,10^4,10^6],'xticklabel',[10^2,10^4,10^6],'fontsize',16); grid on;
    xlabel('Surface area (m^2)','fontsize',16);
    % text(10000,10^-2,'seaward','fontsize',16);
    drawnow;
    
    %plot seasonal profiles of melange characteristics
    % bergAfig = figure; set(gcf,'position',[50 50 1200 400]); 
    figure(sitefig); %subplot(4,2,[7,8]);
    clear pz;
    for k = 1:4
        %site figure
        % figure(bergAfig);
        figure(sitefig); subplot(subp);
        ax2 = gca;
        yyaxis left
        if sum(~isnan(size_plslope(k,:))) > 0
           ps(k) = plot(Hdist,res_A(k,:),'-','color',seas_cmap(k,:),'linewidth',3); hold on; %bergy bit misfits
           %add symbols
           if ismember(MP(j).name,big3)
               plot(Hdist,res_A(k,:),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           else
               plot(Hdist,res_A(k,:),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           end
        end

        %plot the fractional difference between powerlaw-modeled & observed
        %area covered by small icebergs
        plot(Hdist,res_A(k,:),'-','color',seas_cmap(k,:),'linewidth',3); hold on; %bergy bit misfit
        ax2.YLim = [-1,1];

        %plot packing density
        Hmean = nanmean(H_seas(:,k,:),3);
        packmean = 100*nanmean(pack_seas(:,k,:),3);
        if sum(~isnan(Hmean)) > 0
            yyaxis right;
            plot(Hdist(~isnan(Hmean))',packmean(~isnan(Hmean)),'--','color',seas_cmap(k,:),'linewidth',3); hold on;
            %add symbols
           if ismember(MP(j).name,big3)
               plot(Hdist(~isnan(Hmean))',packmean(~isnan(Hmean)),'d','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           else
               plot(Hdist(~isnan(Hmean))',packmean(~isnan(Hmean)),'s','color',seas_cmap(k,:),'linewidth',1,'markersize',5,'markerfacecolor',seas_cmap(k,:)); hold on; %bergy bit misfits
           end
            ax2.YLim = [0,100]; 
        end

        %format the site figure
        ax2.YAxis(1).Color = 'k'; ax2.YAxis(2).Color = 'k'; 
        ax2.YAxis(1).Label.String = 'Bergy bit misfit (%)'; ax2.YAxis(1).Label.FontSize = 16;
        ax2.YAxis(2).Label.String = 'Packing density (%)'; ax2.YAxis(2).Label.FontSize = 16;
        drawnow;

        %scatterplot of bergy bit misfits vs packing density for the site
        figure(sitefig); subplot(subm); set(gca,'box','on');
        Hdist(Hdist==0) = 1;
        if ismember(MP(j).name,big3)
            %plot all data along the profiles
            scatter(-res_A(k,:),packmean,240-18*log(Hdist),'d','MarkerFaceColor','none',...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;
        else
            %plot all data along the profiles
            scatter(-res_A(k,:),packmean,240-18*log(Hdist),'s','MarkerFaceColor','none',...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;
        end
        drawnow;

        %scatterplot of bergy bit misfits vs packing density for all sites
        figure(missfig); set(gca,'box','on');
        if ismember(MP(j).name,big3)
            %plot dummy points for the legend
            small_ref = find(Hdist<1000,1,'first');
            if ~isempty(small_ref) && ~isnan(packmean(small_ref))
                pm(1) = plot(-res_A(k,small_ref),packmean(small_ref),'d','MarkerSize',sqrt(240-18*log(Hdist(small_ref))),'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k','LineWidth',1.5); hold on;
            end
            clear small_ref;
            med_ref = find(Hdist<10000,1,'last');
            if ~isempty(med_ref) && ~isnan(packmean(med_ref))
                pm(2) = plot(-res_A(k,med_ref),packmean(med_ref),'d','MarkerSize',sqrt(240-18*log(Hdist(med_ref))),'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k','LineWidth',1.5); hold on;
            end
            clear med_ref;
            big_ref = find(Hdist<20000,1,'last');
            if ~isempty(big_ref)
                pm(3) = plot(-res_A(k,big_ref),packmean(big_ref),'d','MarkerSize',sqrt(240-18*log(Hdist(big_ref))),'MarkerFaceColor','none',...
                    'MarkerEdgeColor','k','LineWidth',1.5); hold on;
            end
            clear big_ref;

            %plot all data along the profiles
            scatter(-res_A(k,:),packmean,240-18*log(Hdist),'d','MarkerFaceColor','none',...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;

        else
            %plot all data along the profiles
            scatter(-res_A(k,:),packmean,240-18*log(Hdist),'s','MarkerFaceColor','none',...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;
            
        end
        % disp(['Min bergy bit plot symbol size: ',num2str(min(240-18*log(Hdist)))]);
        drawnow; clear packmean Hmean;
    end
    %call out the data for the example map & size distribution plot
    if contains(MP(j).name,'ASG')
        Dsubs = dir([root_dir,MP(j).name,'/',MP(j).name,'*-iceberg-distribution-subsets.csv']);
        for p = 1:length(Dsubs)
            distdate(p,:) = Dsubs(p).name(5:12);
        end
        dateref = strmatch('20110610',distdate);
        %load the data for that date
        D = readtable([root_dir,MP(j).name,'/',Dsubs(dateref).name],"VariableNamingRule","preserve");
        v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
        %first 2 columns are area & bin width, so 3+ are data from points along the centerline
        if ~isnan(inland_idx(p))
            berg_nos = table2array(D(:,(inland_idx(dateref)-1)+2));
            berg_nos(berg_nos==0) = NaN;
            n1 = berg_nos./dv1;
        end
        clear D;
        clear berg_nos Dsubs distdate;
        packing_dated = packing(dateref,inland_idx(dateref)-1);

        %filter to decide if there are enough points to fit lines
        v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
        v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs

        %enter unfiltered data into model (filtering happens there)
        tempfig = figure;
        [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
        n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model
        close(tempfig);

        %calculate misfit for bergy bits
        res = n1_mod-n1; % calculate residuals from E-BC model (data - model)
        res_Adated = sum(v1(v1<=vthresh).*(res(v1<=vthresh).*dv1(v1<=vthresh)))./sum(v1.*(n1.*dv1),"omitnan");
        clear v dv n v1 dv1 n1 dateref res;

        %add example date & site info to the scatterplot
        figure(missfig); 
        %plot the same symbol as all 4 season colors to make a season legend
        for k = 1:4
            scat(k) = scatter(-res_Adated,100*packing_dated,240,'s','MarkerFaceColor',seas_cmap(k,:),...
                'MarkerEdgeColor',seas_cmap(k,:),'LineWidth',1.5); hold on;
        end
        %plot the correct data
        scatter(-res_Adated,100*packing_dated,240,'s','MarkerFaceColor',seas_cmap(2,:),...
            'MarkerEdgeColor',seas_cmap(2,:),'LineWidth',1.5); hold on;
        clear *_dated;
    end

    %format the site figure
    % figure(bergAfig);
    figure(sitefig); subplot(subp);
    % pos = get(gca,'position');
    set(gca,'fontsize',16); set(gca,'xlim',[0,max_xlim]); 
    xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/1000); clear xticks;
    xlabel('Distance from terminus (km)','fontsize',16); 
    yyaxis left
    plot([0,max_xlim],[0,0],'-k','linewidth',1); hold on; %add zero misfit line
    set(gca,'ylim',[-1,1],'ytick',[-1:0.5:1],'yticklabel',[-100:50:100]); %adjust labeling to percents
    subplot(subm);
    set(gca,'fontsize',16); grid on; drawnow;
    % set(gca,'xlim',[-1,1],'xtick',[-1:0.5:1],'xticklabel',[-100:50:100],'ylim',[0,100]);
    set(gca,'ylim',[0,100]);
    plot([0,0],[0,100],'-k','linewidth',1); hold on; %add zero misfit line
    xlabel('Bergy bit misfit (% Area)','fontsize',16); ylabel('Packing density (%)','fontsize',16); 
    %format the bergy bit figure for all sites
    figure(missfig);
    set(gca,'fontsize',16); grid on; drawnow;
    set(gca,'xlim',[-1,1],'xtick',[-1:0.5:1],'xticklabel',[-100:50:100],'ylim',[0,100]);
    xlabel('Bergy bit misfit (% Area)','fontsize',16); ylabel('Packing density (%)','fontsize',16); 
    figure(sitefig);
    leg_pos = get(seas_leg,'Position')
    seas_leg.Position = [leg_pos(1), 0.95 leg_pos(3), leg_pos(4)]; clear leg_pos;
    clear pz ps *_leg pos;

    %compile the buttressing data for plotting
    BM_annual = [BM_annual; squeeze(MP(j).B.butt_Meng(1,1,:))./10^6, squeeze(MP(j).B.butt_Meng(1,2,:))./10^6, squeeze(MP(j).B.butt_Meng(1,3,:))./10^6, squeeze(MP(j).B.butt_Meng(1,4,:))./10^6];
    BM_character = [BM_character; nanmean(MP(j).B.butt_Meng(1,1,:))/10^6, nanmean(MP(j).B.butt_Meng(1,2,:))/10^6, nanmean(MP(j).B.butt_Meng(1,3,:))/10^6, nanmean(MP(j).B.butt_Meng(1,4,:))/10^6];
    BA_annual = [BA_annual; squeeze(MP(j).B.butt_Amundson(1,1,:))./10^6, squeeze(MP(j).B.butt_Amundson(1,2,:))./10^6, squeeze(MP(j).B.butt_Amundson(1,3,:))./10^6, squeeze(MP(j).B.butt_Amundson(1,4,:))./10^6];
    BA_character = [BA_character; nanmean(MP(j).B.butt_Amundson(1,1,:))/10^6, nanmean(MP(j).B.butt_Amundson(1,2,:))/10^6, nanmean(MP(j).B.butt_Amundson(1,3,:))/10^6, nanmean(MP(j).B.butt_Amundson(1,4,:))/10^6];


    %save the data and the figure
    save([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat'],'MP','-v7.3');
    saveas(sitefig,[root_dir,site_abbrev,'/',site_abbrev,'-seasonal-speed-size_',num2str(Hcutoff),'m-Hthreshold_',num2str(vdtmin),'-',num2str(vdtmax),'dt-',vfilter,'-speeds_',sampling,'-profiles.png'],'png'); %save the plots
    exportgraphics(sitefig,[root_dir,site_abbrev,'/',site_abbrev,'-seasonal-speed-size_',num2str(Hcutoff),'m-Hthreshold_',num2str(vdtmin),'-',num2str(vdtmax),'dt-',vfilter,'-speeds_',sampling,'-profiles.tif'],Resolution=600);
    % saveas(bergAfig,[root_dir,MP(j).name,'/',MP(j).name,'-seasonal-melange-properties_profiles.png'],'png'); %save the plots
    disp(['Done with #',num2str(j),': ',site_abbrev]);
    % % figure(sitefig);
    % % disp('Close the site figure to advance'); disp(' ');
    % % uiwait %advance only after figure is closed
    % saveas(missfig,[root_dir,MP(j).name,'/',MP(j).name,'-bergybit_scatterplot.png'],'png'); 
    % exportgraphics(missfig,[root_dir,MP(j).name,'/',MP(j).name,'-bergybit_scatterplot.tif'],Resolution=600);
    % clf(missfig); %clear everything from the bergy bit misfit figure so you can look at data for each site
    close(sitefig); %close(bergAfig);
    drawnow; 

    %refresh
    clear berg_* bergdist* berg_normdist* C centerline* D Dsubs *idx seaward_* inland_* term_* tran_* size_classes Zfilt H_* Havg vel_* vels* v_mean w zdate berg_mo bins bin_no z_* pos pz pv *dist *yrs *mos seas_leg packing;
    clear leg_* ref_* Tdate* xlims Tslopes res* Tres size_plslope bin* column_names dv1 v1 n1* *_seas ax ax2 max_xlim width_prof Vterm dVdxterm subv subz*;
    clear ax1 ax2;

end
figure(missfig); set(gca,'box','on','linewidth',1.5);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
% miss_leg = legend(pm,'terminus','10 km','20 km'); set(miss_leg,'location','southwest');
%spoof a legend for distance from the terminus so it includes both symbol types
rectangle('position',[min(xlims)+0.05*range(xlims), min(ylims)+0.075*range(ylims),0.25*range(xlims),0.3*range(ylims)],'facecolor', ...
    'w','edgecolor','k','linewidth',1.5);
plot(min(xlims)+0.085*range(xlims), min(ylims)+0.18*range(ylims),...
    's','MarkerSize',sqrt(240-18*log(1)),'markeredgecolor','k','markerfacecolor','w','linewidth',1);
plot(min(xlims)+0.085*range(xlims), min(ylims)+0.135*range(ylims),...
    's','MarkerSize',sqrt(240-18*log(10000)),'markeredgecolor','k','markerfacecolor','w','linewidth',1);
plot(min(xlims)+0.085*range(xlims), min(ylims)+0.10*range(ylims),...
    's','MarkerSize',sqrt(240-18*log(20000)),'markeredgecolor','k','markerfacecolor','w','linewidth',1);
plot(min(xlims)+0.14*range(xlims), min(ylims)+0.18*range(ylims),...
    'd','MarkerSize',sqrt(240-18*log(1)),'markeredgecolor','k','markerfacecolor','w','linewidth',1);
plot(min(xlims)+0.14*range(xlims), min(ylims)+0.135*range(ylims),...
    'd','MarkerSize',sqrt(240-18*log(10000)),'markeredgecolor','k','markerfacecolor','w','linewidth',1);
plot(min(xlims)+0.14*range(xlims), min(ylims)+0.10*range(ylims),...
    'd','MarkerSize',sqrt(240-18*log(20000)),'markeredgecolor','k','markerfacecolor','w','linewidth',1);
text(min(xlims)+0.18*range(xlims),min(ylims)+0.185*range(ylims),['terminus'],'fontsize',14);
text(min(xlims)+0.18*range(xlims),min(ylims)+0.14*range(ylims),['10 km'],'fontsize',14);
text(min(xlims)+0.18*range(xlims),min(ylims)+0.1*range(ylims),['20 km'],'fontsize',14);
text(min(xlims)+0.08*range(xlims),min(ylims)+0.22*range(ylims),['D<15 Gt yr^{-1}'],'fontsize',14,'rotation',45);
text(min(xlims)+0.155*range(xlims),min(ylims)+0.22*range(ylims),['D>15 Gt yr^{-1}'],'fontsize',14,'rotation',45);
missleg = legend(scat,season_names); missleg.Location = 'northeast';
saveas(missfig,[root_dir,'GrIS-bergybit_',num2str(zcutoff),'m-zcutoff_scatterplot.png'],'png'); 
exportgraphics(missfig,[root_dir,'GrIS-bergybit_',num2str(zcutoff),'m-zcutoff_scatterplot.tif'],Resolution=600);
clear pm;


%save the GrIS-wide profiles
saveas(summaryfig,[root_dir,'GrIS-melange_thickness-speed_',num2str(zcutoff),'m-zcutoff_profiles.png'],'png'); 
exportgraphics(summaryfig,[root_dir,'GrIS-melange_thickness-speed_',num2str(zcutoff),'m-zcutoff_profiles.tif'],Resolution=600); 


%create histograms of buttressing
buttfig = figure; set(buttfig,'position',[450 50 900 900]); 
for k = [2,3,4,1]
    subplot(2,2,1);
    hM(k) = histogram(BM_annual(~isnan(BM_annual(:,k)),k),'BinEdges',[0.1:0.1:13],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,17],'xlim',[0,13],'box','on'); ylabel('Count','fontsize',16); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.15,0.95*17,'a) annual packing-based buttressing','fontsize',16);
    subplot(2,2,2);
    hA(k) = histogram(BA_annual(~isnan(BA_annual(:,k)),k),'BinEdges',[0.1:0.1:13],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,17],'xlim',[0,13],'box','on'); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.15,0.95*17,'b) annual strain rate-based buttressing','fontsize',16);
    subplot(2,2,3);
    histogram(BM_character(~isnan(BM_character(:,k)),k),'BinEdges',[0.1:0.1:7],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,6],'xlim',[0,7],'box','on'); ylabel('Count','fontsize',16); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.1,0.95*6,'c) characteristic packing-based buttressing','fontsize',16);
    subplot(2,2,4);
    histogram(BA_character(~isnan(BA_character(:,k)),k),'BinEdges',[0.1:0.1:7],'FaceColor',seas_cmap(k,:),'EdgeColor',seas_cmap(k,:),...
        'EdgeAlpha',1,'LineWidth',1); hold on;
    set(gca,'fontsize',16,'ylim',[0,6],'xlim',[0,7],'box','on'); xlabel('Buttressing (x10^6 N/m)','fontsize',16);
    text(0.1,0.95*6,'d) characteristic strain rate-based buttressing','fontsize',16);

    %display medians across all sites for the season
    disp(['Median Meng buttressing estimate for ',char(season_names(k)),' = ',num2str(round(nanmedian(BM_character(:,k),"all"),2)),'x10^6 N/m']);
    disp(['Median Amundson buttressing estimate for ',char(season_names(k)),' = ',num2str(round(nanmedian(BA_character(:,k),"all"),2)),'x10^6 N/m']);
end
subplot(2,2,1); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
subplot(2,2,2); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
subplot(2,2,3); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
subplot(2,2,4); pos = get(gca,'position'); set(gca,'position',[pos(1)-0.03 pos(2) 1.15*pos(3) 1.05*pos(4)]);
disp(['Mean difference between Meng & Amundson buttressing estimates = ',num2str(round(nanmean(BM_annual-BA_annual,"all"),2)),'x10^6 N/m']);
disp(['Mean PERCENT difference between Meng & Amundson buttressing estimates = ',num2str(round(100*nanmean((BM_annual-BA_annual)./nanmean(cat(3,BM_annual,BA_annual),3),"all"),0)),'%']);
%display medians across spring & summer using only sites with data from both seasons
disp('Seasonal Meng & Amundson buttressing omitting sites with spring-only data:')
disp(['Medians for ',char(season_names(2)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)),2),"all"),2)),...
    ' and ',num2str(round(nanmedian(BA_character(:,2),"all"),2)),'x10^6 N/m']);
disp(['Medians for ',char(season_names(3)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)),3),"all"),2)),...
    ' and ',num2str(round(nanmedian([BA_character(~isnan(BA_character(:,3)),3); zeros(sum(~isnan(BM_character(~isnan(BM_character(:,3)),2)))-sum(~isnan(BA_character(:,3))),1)],"all"),2)),'x10^6 N/m']);
disp('Seasonal Meng & Amundson buttressing omitting sites with spring-only data AND without strain rate-based estimates:')
disp(['Medians for ',char(season_names(2)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),2),"all"),2)),...
    ' and ',num2str(round(nanmedian(BA_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),2),"all"),2)),'x10^6 N/m']);
disp(['Medians for ',char(season_names(3)),' = ',num2str(round(nanmedian(BM_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),3),"all"),2)),...
    ' and ',num2str(round(nanmedian(BA_character(~isnan(BM_character(:,3)) & ~isnan(BA_character(:,3)),3),"all"),2)),'x10^6 N/m']);
%save the buttressing plot
saveas(buttfig,[root_dir,'GrIS-melange_buttressing_',num2str(zcutoff),'m-zcutoff_histograms.png'],'png'); 
exportgraphics(buttfig,[root_dir,'GrIS-melange_buttressing_',num2str(zcutoff),'m-zcutoff_histograms.tif'],Resolution=600);


%% plot figures that combine terminus position and buttressing observations
% show the time-stamped buttressing estimates, noting that velocities can
% come from before and after the elevation observations, alongside the
% terminus position timeseries to see if they relate
close all;
% load([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat']);
% for j = 1:length(MP)
%     saved_sites(j,:) = MP(j).name;
% end

%create a composite figure of subplots
summaryfig = figure; set(summaryfig,'position',[-1650 50 900 1200]);
for j = 1:length(plot_locs)
    figure(summaryfig);
    eval(['subB',num2str(geo_ind(j)),'=subplot(',num2str(rows),',',num2str(cols),',',num2str(plot_locs(j)),');']);
end

%create a figure showing near-terminus thickness vs terminus change rate &
%thickness vs buttressing
termfig = figure; set(termfig,'position',[50 50 1800 600]);
sub1 = subplot(1,2,1); sub2 = subplot(1,2,2);

%create a figure that shows terminus change rates around buttressing
%observations for sites with detailed time series
dterm_6mo_fig = figure; set(dterm_6mo_fig,'position',[150 50 900 400]);
dterm_3mo_fig = figure; set(dterm_3mo_fig,'position',[150 150 900 400]);
dterm_1mo_fig = figure; set(dterm_1mo_fig,'position',[150 250 900 400]);
dterm_inst_fig = figure; set(dterm_inst_fig,'position',[150 350 900 800]);

%iterate, adding subplots to the composite & making individual site plots
jams_sr = 0; jams_b = 0; %record the number of jamming events based on compressional strain rates & huge buttressing estimates
dTw = []; dTs = []; Bsp = []; Bsu = []; %dummy matrices for fall to spring terminus change (dTw),spring to fall terminus change (dTs), and spring buttressing (Bsp)
%create dummy vectors to hold seasonal data (individual profiles, not seasonal means)
profs(1).H = NaN(1,50); profs(2).H = NaN(1,50); profs(3).H = NaN(1,50); profs(4).H = NaN(1,50); 
profs(1).p = NaN(1,50); profs(2).p = NaN(1,50); profs(3).p = NaN(1,50); profs(4).p = NaN(1,50); 
profs(1).V = NaN(1,50); profs(2).V = NaN(1,50); profs(3).V = NaN(1,50); profs(4).V = NaN(1,50); 
for l = 1:length(sitenames)
    %identify the index in the structure for the specified sitename
    j = find(strcmp(string(saved_sites),sitenames(l,:)));

    %specify site plotting conventions
    site_abbrev = MP(j).name; disp(site_abbrev);
    if ismember(MP(j).name,big3)
        symbol_shape = 'd';
    else
        symbol_shape = 's';
    end

    %DEM data
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        zdatetime(p) = datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd');
        zyrs(p) = year(zdatetime(p)); zmos(p) = month(zdatetime(p));
        for k = 1:4
            if ismember(zmos(p),seasons(k,:))==1
                zseas(p) = k;
            end
        end
    end
    mean_prof = nanmean(MP(j).Z.transectZavg,2);
    seaward_idx = find(~isnan(mean_prof)==1,1,'first');
    inland_idx = find(~isnan(mean_prof)==1,1,'last');

    %intersect the centerline with each time-stamped melange outline
    for p = 1:length(MP(j).T.date)
        try 
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        catch
            MP(j).T.date{p} = string(str2num(MP(j).T.date{p})-1);
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        end
        Tdate(p) = convert_to_decimaldate(char(MP(j).T.date{p}));
        Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        Tyrs(p) = year(Tdatetime(p)); Tmos(p) = month(Tdatetime(p));
        for k = 1:4
            if ismember(Tmos(p),seasons(k,:))==1
                Tseas(p) = k;
            end
        end
    end

    %compile the data for identifying terminus gaps
    term_decidates = [Tdate, zdate(term_trace == 1)]; term_dates = [Tdatetime, zdatetime(term_trace ==1)]; 
    term_dists = [MP(j).T.termdist, MP(j).Z.termdist(term_trace ==1)];
    term_seas = [Tseas, zseas(term_trace ==1)];
    %sort terminus position observations by date
    [term_decidates_sorted, ic] = sort(term_decidates);
    term_dists_sorted = term_dists(ic); term_seas_sorted = term_seas(ic);
    term_pos_sorted = max(term_dists(term_decidates>min(years)))-term_dists_sorted;

    %compile dated thickness and packing density profiles
    D_dir = [root_dir,MP(j).name,'/'];
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
    end
    [~,Havg,packing,inland_idx,~,~,~,~] = extract_seasonal_binned_DEM_profiles(D_dir,MP(j).name,term_trace,MP(j).Z.termX,MP(j).Z.termY,MP(j).V.X,MP(j).V.Y,ARcomp.best.autoALL,zcutoff,sampling,years,seasons);
    H_profs = NaN(length(zdatetime),max(inland_idx)); pack_profs = NaN(length(zdatetime),max(inland_idx));
    for p = 1:size(Havg,1)
        if isnan(Havg(p,inland_idx(p)-1)) && term_trace(p) == 0
            H_profs(p,1:inland_idx(p)-2) = fliplr(Havg(p,1:inland_idx(p)-2));
            pack_profs(p,1:inland_idx(p)-2) = fliplr(packing(p,1:inland_idx(p)-2));
        else
            H_profs(p,1:inland_idx(p)-1) = fliplr(Havg(p,1:inland_idx(p)-1));
            pack_profs(p,1:inland_idx(p)-1) = fliplr(packing(p,1:inland_idx(p)-1));
        end
    end
    Hdist = MP(j).Z.dist; 

    %compile velocity profiles corresponding to the thickness profiles
    vel_dir = [root_dir,site_abbrev,'/velocities/'];
    if contains(dvterm,'y')
        inland_vel = max(inland_idx); %grab data from the terminus
        ref_adjust = 1; %shift the referencing for the relative coordinate system
    else
        inland_vel = max(inland_idx)-1; %only get speeds from the melange
        ref_adjust = 0;
    end
    vel_profs = create_snapshot_speed_profiles(vel_dir,inland_idx,inland_vel,ref_adjust,vdtmin,vdtmax,zdatetime);
    vdist = MP(j).V.dist; if length(vdist) > size(vel_profs,2); vdist = vdist(1:end-1); end
    dvdist_temp = diff(vdist);
    dVdx = ((diff(vel_profs')'./dvdist_temp)/365);
    %add up apparent jamming events based on compressional strainrates
    if ~isempty(find(dVdx(:,1)<0)); jams_sr = jams_sr + numel(find(dVdx(:,1)<0)); end
    
    %add the profiles for thickness, packing density, and speed to the dummy matrices for each season
    for p = 1:length(zdatetime)
        %isolate seasonal time stamps
        zyrs(p) = year(zdatetime(p)); zmos(p) = month(zdatetime(p));
        [zseas(p),~] = ind2sub(size(seasons),find(seasons == zmos(p)));
    end
    %add to the structures
    Hprofs_dist = [0:1000:49000]; 
    if ~contains(site_abbrev,'ZIM')
        Hdist = [0:nanmedian(diff(MP(j).Z.dist)):(size(H_profs,2)-1)*nanmedian(diff(MP(j).Z.dist))];

        %pad profiles with NaNs on end so they are all the same length
        Hprofs = [H_profs, NaN(size(H_profs,1),length(profs(1).H)-size(H_profs,2))];
        pprofs = [pack_profs, NaN(size(H_profs,1),length(profs(1).H)-size(H_profs,2))];
        % dvprofs = [dVdx, NaN(size(dVdx,1),length(profs(1).H)-size(dVdx,2))];
        for p = 1:length(zdatetime)
            vinterp = interp1(vdist,vel_profs(p,:),Hprofs_dist,'linear');
            if zseas(p) == 1
                profs(1).H = [profs(1).H; Hprofs(p,:)];
                profs(1).p = [profs(1).p; pprofs(p,:)];
                profs(1).V = [profs(1).V; vinterp];
            elseif zseas(p) == 2
                profs(2).H = [profs(2).H; Hprofs(p,:)];
                profs(2).p = [profs(2).p; pprofs(p,:)];
                profs(2).V = [profs(2).V; vinterp];
            elseif zseas(p) == 3
                profs(3).H = [profs(3).H; Hprofs(p,:)];
                profs(3).p = [profs(3).p; pprofs(p,:)];
                profs(3).V = [profs(3).V; vinterp];
            else
                profs(4).H = [profs(4).H; Hprofs(p,:)];
                profs(4).p = [profs(4).p; pprofs(p,:)];
                profs(4).V = [profs(4).V; vinterp];
            end
            clear vinterp;
        end
        clear Hprofs pprofs dvprofs;
    else 
        Hdist = [0:nanmedian(diff(MP(j).Z.dist)):(size(H_profs,2)-1)*nanmedian(diff(MP(j).Z.dist))];
        for p = 1:length(zdatetime)
            %interpolate profiles to 1 km interval (from 1.5 km interval) 
            Hinterp = interp1(Hdist,H_profs(p,:),Hprofs_dist,'linear');
            pinterp = interp1(Hdist,pack_profs(p,:),Hprofs_dist,'linear');
            % dvinterp = interp1((vdist(1:end-1)+vdist(2:end))/2,dVdx(p,:),Hprofs_dist,'linear');
            vinterp = interp1(vdist,vel_profs(p,:),Hprofs_dist,'linear');
            if zseas(p) == 1
                profs(1).H = [profs(1).H; Hinterp];
                profs(1).p = [profs(1).p; pinterp];
                profs(1).V = [profs(1).V; vinterp];
            elseif zseas(p) == 2
                profs(2).H = [profs(2).H; Hinterp];
                profs(2).p = [profs(2).p; pinterp];
                profs(2).V = [profs(2).V; vinterp];
            elseif zseas(p) == 3
                profs(3).H = [profs(3).H; Hinterp];
                profs(3).p = [profs(3).p; pinterp];
                profs(3).V = [profs(3).V; vinterp];
            else
                profs(4).H = [profs(4).H; Hinterp];
                profs(4).p = [profs(4).p; pinterp];
                profs(4).V = [profs(4).V; vinterp];
            end
            clear Hinterp pinterp dvinterp;
        end
    end
    clear Hdist;

    %filter crazy terminus delineations
    dTdt = diff(term_pos_sorted)./diff(term_decidates_sorted);
    %remove duplicates for a given date
    for k = 1:length(dTdt)
        if isinf(dTdt(k))
            term_decidates_sorted(k) = NaN; term_seas_sorted(k) = NaN;
            term_dists_sorted(k+1) = mean(term_dists_sorted(k:k+1)); term_dists_sorted(k) = NaN;
            term_pos_sorted(k+1) = mean(term_pos_sorted(k:k+1)); term_pos_sorted(k) = NaN;
        end
    end
    term_decidates_sorted(isnan(term_pos_sorted)) = []; term_seas_sorted(isnan(term_pos_sorted)) = [];
    term_dists_sorted(isnan(term_pos_sorted)) = []; term_pos_sorted(isnan(term_pos_sorted)) = [];
    clear dTdt;
    dTdt = diff(term_pos_sorted)./diff(term_decidates_sorted);
    term_anoms = find(dTdt > 5*nanmedian(vel_profs(:,1))); %find terminus dates where the next position is way advanced
    bad_flag = []; 
    for k = 1:length(term_anoms)
        if term_anoms(k) == 1
            bad_flag = [bad_flag; term_anoms(k)]; 
        elseif term_anoms(k) > 1 && term_anoms(k) < length(dTdt)
            if dTdt(term_anoms(k)-1) < -5*nanmedian(vel_profs(:,1)) %if the previous change was also a big retreat, its a dip to remove
                bad_flag = [bad_flag; term_anoms(k)];
            elseif dTdt(term_anoms(k)+1) < -5*nanmedian(vel_profs(:,1)) %if the next change was also a big retreat, its a jump to remove
                bad_flag = [bad_flag; term_anoms(k)+1];
            end
        else
            if dTdt(term_anoms(k)-1) < -5*nanmedian(vel_profs(:,1)) %if the previous change was also a big retreat, its a dip to remove
                bad_flag = [bad_flag; term_anoms(k)];
            end
        end
    end
    term_decidates_sorted(bad_flag) = []; term_seas_sorted(bad_flag) = [];
    term_dists_sorted(bad_flag) = []; term_pos_sorted(bad_flag) = []; clear term_anoms bad_flag dTdt;

    %calculate buttressing for each observation
    dTdt_6mo_pre = []; dTdt_6mo_post = []; dTdt_3mo_pre = []; dTdt_3mo_post = []; dTdt_1mo_pre = []; dTdt_1mo_post = []; dTdt_pre = []; dTdt_post = []; %dummy matrices for for terminus change rates before & after thickness observations
    seas_6mo_pre = []; seas_6mo_post = [];
    HT = []; BT = []; dBT_pre = []; dBT_post = []; %dummy matrices for thickness & buttressing estimates
    for p = 1:size(Havg,1)
        HT = [HT; H_profs(p,1)]; 
        press = 0.5*rho_i*(1-(rho_i/rho_sw))*9.81*H_profs(p,1);
        butt_Meng(p) = press*pack_profs(p,1)*H_profs(p,1); BT = [BT; butt_Meng(p)];
        butt_Amundson(p) = (-2*(H_profs(p,1)*press*dVdx(p,1))/((dVdx(p,1)/0.3)+dVdx(p,1)))+press*H_profs(p,1);
        clear press dvdist_temp;

        %calculate terminus position change for a few different time periods
        %relative to each buttressing observation: +/-6 mo, +/-3 mo, +/-1 mo, &
        %between buttressing observations
        %NOTE: minus sign is to make advance positive
        if term_trace(p) == 1
            %6 mo
            before_idx = find(term_decidates_sorted>zdate(p)-0.5 & term_decidates_sorted<zdate(p));
            after_idx = find(term_decidates_sorted>zdate(p) & term_decidates_sorted<zdate(p)+0.5);
            if ~isempty(before_idx)
                dTdt_6mo_pre = [dTdt_6mo_pre; -((MP(j).Z.termdist(p)-term_dists_sorted(before_idx(1)))/abs(zdate(p)-term_decidates_sorted(before_idx(1))))/365]; 
            else
                dTdt_6mo_pre = [dTdt_6mo_pre; NaN];
            end
            if ~isempty(after_idx)
                dTdt_6mo_post = [dTdt_6mo_post; -((term_dists_sorted(after_idx(end))-MP(j).Z.termdist(p))/abs(term_decidates_sorted(after_idx(end))-zdate(p)))/365]; 
            else
                dTdt_6mo_post = [dTdt_6mo_post; NaN];
            end
            seas_6mo_pre = [seas_6mo_pre; zseas(p)]; seas_6mo_post = [seas_6mo_post; zseas(p)]; 
            clear *_idx;

            %3 mo
            before_idx = find(term_decidates_sorted>zdate(p)-0.25 & term_decidates_sorted<zdate(p));
            after_idx = find(term_decidates_sorted>zdate(p) & term_decidates_sorted<zdate(p)+0.25);
            if ~isempty(before_idx)
                dTdt_3mo_pre = [dTdt_3mo_pre; -((MP(j).Z.termdist(p)-term_dists_sorted(before_idx(1)))/abs(zdate(p)-term_decidates_sorted(before_idx(1))))/365]; 
            else
                dTdt_3mo_pre = [dTdt_3mo_pre; NaN];
            end
            if ~isempty(after_idx)
                dTdt_3mo_post = [dTdt_3mo_post; -((term_dists_sorted(after_idx(end))-MP(j).Z.termdist(p))/abs(term_decidates_sorted(after_idx(end))-zdate(p)))/365]; 
            else
                dTdt_3mo_post = [dTdt_3mo_post; NaN];
            end
            clear *_idx;

            %1 mo
            before_idx = find(term_decidates_sorted>zdate(p)-0.0833 & term_decidates_sorted<zdate(p));
            after_idx = find(term_decidates_sorted>zdate(p) & term_decidates_sorted<zdate(p)+0.0833);
            if ~isempty(before_idx)
                dTdt_1mo_pre = [dTdt_1mo_pre; -((MP(j).Z.termdist(p)-term_dists_sorted(before_idx(1)))/abs(zdate(p)-term_decidates_sorted(before_idx(1))))/365]; 
            else
                dTdt_1mo_pre = [dTdt_1mo_pre; NaN];
            end
            if ~isempty(after_idx)
                dTdt_1mo_post = [dTdt_1mo_post; -((term_dists_sorted(after_idx(end))-MP(j).Z.termdist(p))/abs(term_decidates_sorted(after_idx(end))-zdate(p)))/365]; 
            else
                dTdt_1mo_post = [dTdt_1mo_post; NaN];
            end
            clear *_idx;

            %between melange observations
            if p > 1
                dTdt_pre = [dTdt_pre; -((MP(j).Z.termdist(p)-MP(j).Z.termdist(p-1))/abs(zdate(p)-zdate(p-1)))/365]; 
            else
                dTdt_pre = [dTdt_pre; NaN];
            end
            if p < size(Havg,1)
                dTdt_post = [dTdt_post; -((MP(j).Z.termdist(p+1)-MP(j).Z.termdist(p))/abs(zdate(p+1)-zdate(p)))/365]; 
            else
                dTdt_post = [dTdt_post; NaN];
            end

        else
            dTdt_6mo_pre = [dTdt_6mo_pre; NaN]; dTdt_6mo_post = [dTdt_6mo_post; NaN];
            seas_6mo_pre = [seas_6mo_pre; NaN]; seas_6mo_post = [seas_6mo_post; NaN]; 
            dTdt_3mo_pre = [dTdt_3mo_pre; NaN]; dTdt_3mo_post = [dTdt_3mo_post; NaN];
            dTdt_1mo_pre = [dTdt_1mo_pre; NaN]; dTdt_1mo_post = [dTdt_1mo_post; NaN];
            dTdt_pre = [dTdt_pre; NaN]; dTdt_post = [dTdt_post; NaN];
        end

    end
    dBT_pre = [NaN; BT(2:end)-BT(1:end-1)]; dBT_post = [BT(2:end)-BT(1:end-1); NaN];

    %plot the buttressing & terminus change data for ASG & ZIM
    if contains(site_abbrev,'ASG') || contains(site_abbrev,'ZIM')
        %specify outlier percentile thresholds & marker sizes
        outlier_percs = [1,99]; good_size = 50; bad_size = 25;

        %specify plot limits
        xlims = [0, 2]; %buttressing range (10^6 N/m)
        if contains(site_abbrev,'ASG')
            ylims = [-40, 40]; yticks = [-40:20:40]; %terminus change rate range (m/d)
        else
            ylims = [-240, 240]; yticks = [-240:80:240]; %terminus change rate range (m/d)
        end

        %plot terminus change +/-6 mo from buttressing obs v buttressing
        figure(dterm_6mo_fig);
        %before
        if contains(site_abbrev,'ASG'); subplot(1,2,1); else; subplot(1,2,2); end
        for k = 1:4
            sp(k) = scatter(BT(1)/10^6,dTdt_6mo_pre(1),good_size,seas_cmap(k,:),'filled',symbol_shape); hold on;
        end
        % scatter(BT/10^6,dTdt_6mo_pre,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        scatter(BT/10^6,dTdt_6mo_pre,good_size,seas_cmap(zseas,:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        dTdtf = isoutlier(dTdt_6mo_pre,"percentiles",outlier_percs);
        % scatter(BT(~isnan(dTdt_6mo_pre) & ~dTdtf)/10^6,dTdt_6mo_pre(~isnan(dTdt_6mo_pre) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_6mo_pre) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        [f,gof] = fit(BT(~isnan(dTdt_6mo_pre) & ~dTdtf)/10^6,dTdt_6mo_pre(~isnan(dTdt_6mo_pre) & ~dTdtf),'poly1');
        disp(['Change in terminus change rate with buttressing = ',num2str(f.p1),' m/d per 10^6 N/m']);
        disp(['R^2 btw term change rate 6 mo before & buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        if gof.rsquare > 0.5; plot([min(BT) max(BT)]/10^6,feval(f,[min(BT) max(BT)]/10^6),'--k'); end
        % set(gca,'fontsize',16,'xlim',xlims,'xticklabel',[],'ylim',ylims,'ytick',yticks,'box','on'); grid on;
        if contains(site_abbrev,'ASG')
            % text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'a) 6 months prior','fontsize',14);
            text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'a) Nunatakassaap','fontsize',16);
            ylabel('Terminus change rate (m/d)','fontsize',16);
            sp_leg = legend(sp,season_names,'location','northoutside','orientation','horizontal','fontsize',16); sp_leg.Position = [0.43 0.93 0.2 0.05];
            pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.1*pos(3) 0.75]);
        else
            % text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['b) 6 months prior'],'fontsize',14);
            text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['b) Zachari',char(230),' Isstr',char(248),'m'],'fontsize',16);
            pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.1*pos(3) 0.75]);
        end
        clear dTdtf f gof;
        % %after 
        % if contains(site_abbrev,'ASG'); subplot(2,2,3); else; subplot(2,2,4); end
        % scatter(BT/10^6,dTdt_6mo_post,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        % dTdtf = isoutlier(dTdt_6mo_post,"percentiles",outlier_percs);
        % scatter(BT(~isnan(dTdt_6mo_post) & ~dTdtf)/10^6,dTdt_6mo_post(~isnan(dTdt_6mo_post) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_6mo_post) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        % [f,gof] = fit(BT(~isnan(dTdt_6mo_post) & ~dTdtf)/10^6,dTdt_6mo_post(~isnan(dTdt_6mo_post) & ~dTdtf),'poly1');
        % disp(['R^2 btw term change rate 6 mo after & buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        % if gof.rsquare > 0.5; plot([min(BT) max(BT)]/10^6,feval(f,[min(BT) max(BT)]/10^6),'--k'); end
        set(gca,'fontsize',16,'xlim',xlims,'ylim',ylims,'ytick',yticks,'box','on'); xlabel('Buttressing (10^6 N/m)','fontsize',16); grid on; clear dTdtf f gof;
        % if contains(site_abbrev,'ASG')
        %     text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'c) 6 months after','fontsize',14);
        %     ylabel('Terminus change rate (m/d)','fontsize',16); 
        % else
        %     text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['d) 6 months after'],'fontsize',14);
        %     annotation(dterm_6mo_fig,'textbox',[0.24,0.9,0.2,0.05],'string','Nunatakassaap','edgecolor','none','fontsize',16);
        %     annotation(dterm_6mo_fig,'textbox',[0.69,0.9,0.2,0.05],'string',['Zachari',char(230),' Isstr',char(248),'m'],'edgecolor','none','fontsize',16);
        % end
        % clear dTdtf f gof;
        % pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.1*pos(3) 1.1*pos(4)]);

        %for 6 month dataset, calculate fits & percent retreating
        %observations for just the spring
        bigB_ind = find(BT>=0.5e6); smallB_ind = find(BT<0.5e6); 
        bigB_dTdt = dTdt_6mo_pre(bigB_ind); smallB_dTdt = dTdt_6mo_pre(smallB_ind); 
        disp(['Fraction with net advance for low buttressing = ',num2str(numel(smallB_dTdt(smallB_dTdt>0))/numel(smallB_dTdt))]);
        disp(['Fraction with net advance for high buttressing = ',num2str(numel(bigB_dTdt(bigB_dTdt>0))/numel(bigB_dTdt))]);
        clear bigB* smallB*;
        BT_sub = BT(zseas(~isnan(dTdt_6mo_pre))==2); dTdt_6mo_pre_sub = dTdt_6mo_pre(zseas(~isnan(dTdt_6mo_pre))==2);
        dTdtf = isoutlier(dTdt_6mo_pre_sub,"percentiles",outlier_percs); 
        [f,gof] = fit(BT_sub(~isnan(dTdt_6mo_pre_sub) & ~dTdtf)/10^6,dTdt_6mo_pre_sub(~isnan(dTdt_6mo_pre_sub) & ~dTdtf),'poly1');
        disp(['Change in terminus change rate with SPRING buttressing = ',num2str(f.p1),' m/d per 10^6 N/m']);
        disp(['R^2 btw term change rate 6 mo before & SPRING buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        bigB_ind = find(BT(zseas == 2)>=0.5e6); smallB_ind = find(BT(zseas == 2)<0.5e6); 
        bigB_dTdt = dTdt_6mo_pre(zseas == 2); bigB_dTdt = bigB_dTdt(bigB_ind); 
        smallB_dTdt = dTdt_6mo_pre(zseas == 2); smallB_dTdt = smallB_dTdt(smallB_ind); 
        disp(['Fraction with net advance for low SPRING buttressing = ',num2str(numel(smallB_dTdt(smallB_dTdt>0))/numel(smallB_dTdt))]);
        disp(['Fraction with net advance for high SPRING buttressing = ',num2str(numel(bigB_dTdt(bigB_dTdt>0))/numel(bigB_dTdt))]);
        clear BT_sub dTdt_*_sub bigB* smallB*;


        %plot terminus change +/-3 mo from buttressing obs v buttressing
        figure(dterm_3mo_fig);
        %before
        if contains(site_abbrev,'ASG'); subplot(1,2,1); else; subplot(1,2,2); end
        % scatter(BT/10^6,dTdt_3mo_pre,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        scatter(BT/10^6,dTdt_3mo_pre,good_size,seas_cmap(zseas,:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        dTdtf = isoutlier(dTdt_3mo_pre,"percentiles",outlier_percs);
        % scatter(BT(~isnan(dTdt_3mo_pre) & ~dTdtf)/10^6,dTdt_3mo_pre(~isnan(dTdt_3mo_pre) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_3mo_pre) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        [f,gof] = fit(BT(~isnan(dTdt_3mo_pre) & ~dTdtf)/10^6,dTdt_3mo_pre(~isnan(dTdt_3mo_pre) & ~dTdtf),'poly1');
        disp(['Change in terminus change rate with buttressing = ',num2str(f.p1),' m/d per 10^6 N/m']);
        disp(['R^2 btw term change rate 3 mo before & buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        if gof.rsquare > 0.5; plot([min(BT) max(BT)]/10^6,feval(f,[min(BT) max(BT)]/10^6),'--k'); end
        % set(gca,'fontsize',16,'xlim',xlims,'xticklabel',[],'ylim',ylims,'ytick',yticks,'box','on'); grid on; clear dTdtf f gof;
        if contains(site_abbrev,'ASG')
            % text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'a) 3 months prior','fontsize',14);
            text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'a) Nunatakassaap','fontsize',16);
            ylabel('Terminus change rate (m/d)','fontsize',16); 
        else
            % text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['b) 3 months prior'],'fontsize',14);
            text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['b) Zachari',char(230),' Isstr',char(248),'m'],'fontsize',16);
        end
        clear dTdtf f gof;
        pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.02 1.1*pos(3) pos(4)]);
        % %after 
        % if contains(site_abbrev,'ASG'); subplot(2,2,3); else; subplot(2,2,4); end
        % scatter(BT/10^6,dTdt_3mo_post,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        % dTdtf = isoutlier(dTdt_3mo_post,"percentiles",outlier_percs);
        % scatter(BT(~isnan(dTdt_3mo_post) & ~dTdtf)/10^6,dTdt_3mo_post(~isnan(dTdt_3mo_post) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_3mo_post) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        % [f,gof] = fit(BT(~isnan(dTdt_3mo_post) & ~dTdtf)/10^6,dTdt_3mo_post(~isnan(dTdt_3mo_post) & ~dTdtf),'poly1');
        % disp(['R^2 btw term change rate 3 mo after & buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        % if gof.rsquare > 0.5; plot([min(BT) max(BT)]/10^6,feval(f,[min(BT) max(BT)]/10^6),'--k'); end
        set(gca,'fontsize',16,'xlim',xlims,'ylim',ylims,'ytick',yticks,'box','on'); xlabel('Buttressing (10^6 N/m)','fontsize',16); grid on; clear dTdtf f gof;
        % if contains(site_abbrev,'ASG')
        %     text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'c) 3 months after','fontsize',14);
        %     ylabel('Terminus change rate (m/d)','fontsize',16); 
        % else
        %     text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['d) 3 months after'],'fontsize',14);
        %     annotation(dterm_3mo_fig,'textbox',[0.24,0.9,0.2,0.05],'string','Nunatakassaap','edgecolor','none','fontsize',16);
        %     annotation(dterm_3mo_fig,'textbox',[0.69,0.9,0.2,0.05],'string',['Zachari',char(230),' Isstr',char(248),'m'],'edgecolor','none','fontsize',16);
        % end
        % clear dTdtf f gof;
        % pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.1*pos(3) 1.1*pos(4)]);

        %plot terminus change +/-1 mo from buttressing obs v buttressing
        figure(dterm_1mo_fig);
        %before
        if contains(site_abbrev,'ASG'); subplot(1,2,1); else; subplot(1,2,2); end
        % scatter(BT/10^6,dTdt_1mo_pre,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        scatter(BT/10^6,dTdt_1mo_pre,good_size,seas_cmap(zseas,:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        dTdtf = isoutlier(dTdt_1mo_pre,"percentiles",outlier_percs);
        % scatter(BT(~isnan(dTdt_1mo_pre) & ~dTdtf)/10^6,dTdt_1mo_pre(~isnan(dTdt_1mo_pre) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_1mo_pre) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        [f,gof] = fit(BT(~isnan(dTdt_1mo_pre) & ~dTdtf)/10^6,dTdt_1mo_pre(~isnan(dTdt_1mo_pre) & ~dTdtf),'poly1');
        disp(['Change in terminus change rate with buttressing = ',num2str(f.p1),' m/d per 10^6 N/m']);
        disp(['R^2 btw term change rate 1 mo before & buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        if gof.rsquare > 0.5; plot([min(BT) max(BT)]/10^6,feval(f,[min(BT) max(BT)]/10^6),'--k'); end
        % set(gca,'fontsize',16,'xlim',xlims,'xticklabel',[],'ylim',ylims,'ytick',yticks,'box','on'); grid on; clear dTdtf f gof;
        if contains(site_abbrev,'ASG')
            % text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'a) 1 month prior','fontsize',14);
            text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'a) Nunatakassaap','fontsize',16);
            ylabel('Terminus change rate (m/d)','fontsize',16); 
        else
            % text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['b)
            % 1 month prior'],'fontsize',14);e
            text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['b) Zachari',char(230),' Isstr',char(248),'m'],'fontsize',16);
        end
        clear dTdtf f gof;
        pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.02 1.1*pos(3) pos(4)]);
        % %after 
        % if contains(site_abbrev,'ASG'); subplot(2,2,3); else; subplot(2,2,4); end
        % scatter(BT/10^6,dTdt_1mo_post,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        % dTdtf = isoutlier(dTdt_1mo_post,"percentiles",outlier_percs);
        % scatter(BT(~isnan(dTdt_1mo_post) & ~dTdtf)/10^6,dTdt_1mo_post(~isnan(dTdt_1mo_post) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_1mo_post) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        % [f,gof] = fit(BT(~isnan(dTdt_1mo_post) & ~dTdtf)/10^6,dTdt_1mo_post(~isnan(dTdt_1mo_post) & ~dTdtf),'poly1');
        % disp(['R^2 btw term change rate 1 mo after & buttressing = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        % if gof.rsquare > 0.5; plot([min(BT) max(BT)]/10^6,feval(f,[min(BT) max(BT)]/10^6),'--k'); end
        set(gca,'fontsize',16,'xlim',xlims,'ylim',ylims,'ytick',yticks,'box','on'); xlabel('Buttressing (10^6 N/m)','fontsize',16); grid on; clear dTdtf f gof;
        % if contains(site_abbrev,'ASG')
        %     text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),'c) 1 month after','fontsize',14);
        %     ylabel('Terminus change rate (m/d)','fontsize',16); 
        % else
        %     text(min(xlims)+0.025*range(xlims),min(ylims)+0.95*range(ylims),['d) 1 month after'],'fontsize',14);
        %     annotation(dterm_1mo_fig,'textbox',[0.24,0.9,0.2,0.05],'string','Nunatakassaap','edgecolor','none','fontsize',16);
        %     annotation(dterm_1mo_fig,'textbox',[0.69,0.9,0.2,0.05],'string',['Zachari',char(230),' Isstr',char(248),'m'],'edgecolor','none','fontsize',16);
        % end
        % clear dTdtf f gof;
        % pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.1*pos(3) 1.1*pos(4)]);

        %plot terminus change between buttressing obs v buttressing change
        figure(dterm_inst_fig);
        %before
        if contains(site_abbrev,'ASG'); subplot(2,2,1); else; subplot(2,2,2); end
        scatter(dBT_pre/10^6,dTdt_pre,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        dTdtf = isoutlier(dTdt_pre,"percentiles",outlier_percs);
        scatter(dBT_pre(~isnan(dTdt_pre) & ~dTdtf)/10^6,dTdt_pre(~isnan(dTdt_pre) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_pre) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        [f,gof] = fit(dBT_pre(~isnan(dTdt_pre) & ~dTdtf & ~isnan(dBT_pre))/10^6,dTdt_pre(~isnan(dTdt_pre) & ~dTdtf & ~isnan(dBT_pre)),'poly1');
        disp(['R^2 btw term change rate immediately before & buttressing change = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        if gof.rsquare > 0.5; plot([min(dBT_pre) max(dBT_pre)]/10^6,feval(f,[min(dBT_pre) max(dBT_pre)]/10^6),'--k'); end
        set(gca,'fontsize',16,'xlim',[-max(xlims),max(xlims)],'xticklabel',[],'ylim',ylims,'box','on'); grid on; clear dTdtf f gof;
        if contains(site_abbrev,'ASG')
            text(-max(xlims)+0.05*range(xlims),min(ylims)+0.95*range(ylims),'a) since last DEM','fontsize',14);
            ylabel('Terminus change rate (m/d)','fontsize',16); 
        else
            text(-max(xlims)+0.05*range(xlims),min(ylims)+0.95*range(ylims),['b) since last DEM'],'fontsize',14);
        end
        clear dTdtf f gof;
        pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.05 1.1*pos(3) 1.1*pos(4)]);
        %after 
        if contains(site_abbrev,'ASG'); subplot(2,2,3); else; subplot(2,2,4); end
        scatter(dBT_post/10^6,dTdt_post,bad_size,seas_cmap(zseas,:),'filled',symbol_shape); hold on;
        dTdtf = isoutlier(dTdt_post,"percentiles",outlier_percs);
        scatter(dBT_post(~isnan(dTdt_post) & ~dTdtf)/10^6,dTdt_post(~isnan(dTdt_post) & ~dTdtf),good_size,seas_cmap(zseas(~isnan(dTdt_post) & ~dTdtf),:),'filled',symbol_shape,'markeredgecolor','k'); hold on; 
        [f,gof] = fit(dBT_post(~isnan(dTdt_post) & ~dTdtf & ~isnan(dBT_post))/10^6,dTdt_post(~isnan(dTdt_post) & ~dTdtf & ~isnan(dBT_post)),'poly1');
        disp(['R^2 btw term change rate immediately after & buttressing change = ',num2str(gof.rsquare),' (',num2str(sum(dTdtf)),' outliers)']);
        if gof.rsquare > 0.5; plot([min(dBT_post) max(dBT_post)]/10^6,feval(f,[min(dBT_post) max(dBT_post)]/10^6),'--k'); end
        set(gca,'fontsize',16,'xlim',[-max(xlims),max(xlims)],'ylim',ylims,'box','on'); xlabel('Buttressing change (10^6 N/m)','fontsize',16); grid on; clear dTdtf f gof;
        if contains(site_abbrev,'ASG')
            text(-max(xlims)+0.05*range(xlims),min(ylims)+0.95*range(ylims),'c) until next DEM','fontsize',14);
            ylabel('Terminus change rate (m/d)','fontsize',16); 
        else
            text(-max(xlims)+0.05*range(xlims),min(ylims)+0.95*range(ylims),['d) until next DEM'],'fontsize',14);
            annotation(dterm_inst_fig,'textbox',[0.24,0.9,0.2,0.05],'string','Nunatakassaap','edgecolor','none','fontsize',16);
            annotation(dterm_inst_fig,'textbox',[0.69,0.9,0.2,0.05],'string',['Zachari',char(230),' Isstr',char(248),'m'],'edgecolor','none','fontsize',16);
        end
        clear dTdtf f gof;
        pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.1*pos(3) 1.1*pos(4)]);

        %save the figures
        saveas(dterm_6mo_fig,[root_dir,'GrIS-melange_6mo-termchange-buttressing_',num2str(zcutoff),'m-zcutoff_scatterplots.png'],'png');
        exportgraphics(dterm_6mo_fig,[root_dir,'GrIS-melange_6mo-termchange-buttressing_',num2str(zcutoff),'m-zcutoff_scatterplots.tif'],Resolution=600);
        saveas(dterm_3mo_fig,[root_dir,'GrIS-melange_3mo-termchange-buttressing_',num2str(zcutoff),'m-zcutoff_scatterplots.png'],'png');
        exportgraphics(dterm_3mo_fig,[root_dir,'GrIS-melange_3mo-termchange-buttressing_',num2str(zcutoff),'m-zcutoff_scatterplots.tif'],Resolution=600);
        saveas(dterm_1mo_fig,[root_dir,'GrIS-melange_1mo-termchange-buttressing_',num2str(zcutoff),'m-zcutoff_scatterplots.png'],'png');
        exportgraphics(dterm_1mo_fig,[root_dir,'GrIS-melange_1mo-termchange-buttressing_',num2str(zcutoff),'m-zcutoff_scatterplots.tif'],Resolution=600);
        saveas(dterm_inst_fig,[root_dir,'GrIS-melange_termchange-buttressingchange_',num2str(zcutoff),'m-zcutoff_scatterplots.png'],'png');
        exportgraphics(dterm_inst_fig,[root_dir,'GrIS-melange_termchange-buttressingchange_',num2str(zcutoff),'m-zcutoff_scatterplots.tif'],Resolution=600);
        % clear dtermfig;
    end

    %add up apparent jamming events based on buttressing
    if ~isempty(find(butt_Meng>1e7)); jams_b = jams_b + numel(find(butt_Meng>1e7)); end

    %compile thickness and buttressing data for each year
    Bseas_med = []; Bseas_mad = []; Hseas_med = []; Hseas_mad = []; Tseas_med = []; Tseas_mad = []; %dummy matrices for buttressing, thickness, and terminus position
    Iseas = []; Sseas = []; %dummy matrices for the site index (I), season (S)
    for p = 1:length(years)
        Iseas = [Iseas; repmat(j,1,4)]; Sseas = [Sseas; 1 2 3 4]; 

        %compile buttressing and associated terminus position and thickness data
        yr_idx = find(zyrs == years(p));
        Bmed = NaN(1,4); Bmad = NaN(1,4); Hmed = NaN(1,4); Hmad = NaN(1,4); Tmed = NaN(1,4); Tmad = NaN(1,4); 
        if ~isempty(yr_idx)
            mos_yr = zmos(yr_idx);
            for k = 1:4
                %data for each season (if it exists, NaN if not)
                if sum(ismember(mos_yr,seasons(k,:))) > 0
                    Bmed(k) = median(butt_Meng(yr_idx(ismember(mos_yr,seasons(k,:))==1)),"all","omitmissing");
                    Bmad(k) = mad(butt_Meng(yr_idx(ismember(mos_yr,seasons(k,:))==1)),1,"all");
                    Hmed(k) = median(HT(yr_idx(ismember(mos_yr,seasons(k,:))==1)),"all","omitmissing");
                    Hmad(k) = mad(HT(yr_idx(ismember(mos_yr,seasons(k,:))==1)),1,"all");
                    Tmed(k) = max(term_dists(term_decidates>min(years))) - median(MP(j).Z.termdist(yr_idx(ismember(mos_yr,seasons(k,:))==1)),"all","omitmissing");
                    Tmad(k) = mad(MP(j).Z.termdist(yr_idx(ismember(mos_yr,seasons(k,:))==1)),1,"all");
                end
            end
        end
        Bseas_med = [Bseas_med; Bmed];
        Bseas_mad = [Bseas_mad; Bmad];
        Hseas_med = [Hseas_med; Hmed];
        Hseas_mad = [Hseas_mad; Hmad];
        Tseas_med = [Tseas_med; Tmed];
        Tseas_mad = [Tseas_mad; Tmad];
        Bsp = [Bsp; Bmed(2)]; Bsu = [Bsu; Bmed(3)];
        clear yr_idx mos_yr;

        %compile terminus position change for the preceding winter and
        %subsequent summer: uses POSITION relative to the terminus, which
        %is the max(distance) - distance, not the distance like in the monthly calculations above
        prevyr_idx = find(floor(term_decidates_sorted) == years(p)-1);
        curryr_idx = find(floor(term_decidates_sorted) == years(p)); currsp_idx = find(term_seas_sorted(curryr_idx)==2);
        [prevT_min,prevT_min_idx] = min(term_pos_sorted(prevyr_idx)); % minimum (most retreated) terminus position for previous year
        % [currT_max,currT_max_idx] = max(term_pos_sorted(curryr_idx(currsp_idx))); % maximum (most advanced) SPRING terminus position for current year
        currT_max = Tseas_med(p,2);
        currT_min = min(term_pos_sorted(curryr_idx)); % minimum (most retreated) terminus position for current year
        dTw_temp = currT_max-prevT_min; if isempty(dTw_temp); dTw_temp = NaN; end
        dTs_temp = currT_min - currT_max; if isempty(dTs_temp); dTs_temp = NaN; end
        dTw = [dTw; dTw_temp]; dTs = [dTs; dTs_temp]; %advance > 0, retreat < 0
        clear *_idx *_min *_max dT*_temp;
    end
    % if j == 15
    %     disp(['Winter terminus change = ',num2str(dTw(end-length(years)+1:end)')]);
    %     disp(['Summer terminus change = ',num2str(dTs(end-length(years)+1:end)')]);
    % end

    %plot the thickness vs terminus change & buttressing
    figure(termfig);
    subplot(sub1); semilogy(H_profs(:,1),butt_Meng,'sk'); hold on;
    for p = 1:size(Havg,1)
        if term_trace(p) == 1
            before_idx = find(term_decidates_sorted<zdate(p),1,'last');
            after_idx = find(term_decidates_sorted>zdate(p),1,'first');

            %dummy plot for legend
            subplot(sub2);
            if j == 1 && p == 1
                pb(1) = semilogx(butt_Meng(p),(-(MP(j).Z.termdist(p)-term_dists_sorted(before_idx))/abs(zdate(p)-term_decidates_sorted(before_idx)))/365,...
                    '^k','markerfacecolor','w'); hold on; %minus sign to make retreat negative
                pb(2) = semilogx(butt_Meng(p),((MP(j).Z.termdist(p)-term_dists_sorted(after_idx))/abs(zdate(p)-term_decidates_sorted(after_idx)))/365,...
                    'vk','markerfacecolor','w'); hold on; %later date minus earlier date so no negative needed
            end

            %plot real data
            subplot(sub2);
            semilogx(butt_Meng(p),(-(MP(j).Z.termdist(p)-term_dists_sorted(before_idx))/abs(zdate(p)-term_decidates_sorted(before_idx)))/365,...
                '^k','markerfacecolor',seas_cmap(zseas(p),:)); hold on; %minus sign to make retreat negative
            semilogx(butt_Meng(p),((MP(j).Z.termdist(p)-term_dists_sorted(after_idx))/abs(zdate(p)-term_decidates_sorted(after_idx)))/365,...
                'vk','markerfacecolor',seas_cmap(zseas(p),:)); hold on; %later date minus earlier date so no negative needed
            clear before_idx after_idx;
        end
    end
    pb_leg = legend(pb,'before','after');
    drawnow;

    %plot the terminus & buttressing timeseries as a standalone site figure
    ts_fig = figure; set(ts_fig,'position',[50 50 1200 600]);
    for k = 1:4
        seas_idx = find(zseas==k);
        for p = 1:length(seas_idx)
            bpM = fill([zdate(seas_idx(p))-0.02,zdate(seas_idx(p))+0.02,zdate(seas_idx(p))+0.02,zdate(seas_idx(p))-0.02,zdate(seas_idx(p))-0.02],...
                [0,0,butt_Meng(seas_idx(p)),butt_Meng(seas_idx(p)),0]./10^6,'w',...
                'edgecolor',seas_cmap(k,:),'facealpha',0,'linewidth',1.5); hold on;
        end
    end
    for k = 1:4
        seas_idx = find(zseas==k);
        for p = 1:length(seas_idx)
            bpA = fill([zdate(seas_idx(p))-0.02,zdate(seas_idx(p))+0.02,zdate(seas_idx(p))+0.02,zdate(seas_idx(p))-0.02,zdate(seas_idx(p))-0.02],...
                [0,0,butt_Amundson(seas_idx(p)),butt_Amundson(seas_idx(p)),0]./10^6,seas_cmap(k,:),...
                'edgecolor','k','facealpha',0.5,'linewidth',0.5); hold on;
        end
    end
    set(gca,'xlim',[min(years)-1 max(years)],'fontsize',16); grid on;
    xlabel('Year','fontsize',16); ylabel('Buttressing (10^6 N/m)','fontsize',16); 

    %overlay the terminus position timeseries
    yyaxis right; axr = gca;
    plot(term_decidates_sorted,term_pos_sorted,'-k','linewidth',1,'markersize',10); hold on;
    for p = 1:length(term_pos_sorted)
        %add dummy plot for the legend
        if p ==1
            for k = 1:4
            pterm(k) = plot(term_decidates_sorted(p),term_pos_sorted(p),symbol_shape,...
            'markeredgecolor',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'linewidth',1,'markersize',10); hold on;
            end
        end

        %plot the data
        plot(term_decidates_sorted(p),term_pos_sorted(p),symbol_shape,...
            'color',seas_cmap(term_seas_sorted(p),:),'markerfacecolor',seas_cmap(term_seas_sorted(p),:),'linewidth',1,'markersize',10); hold on;
    end
    % for p = 1:length(MP(j).Z.date)
    %     plot(zdate(p),max(term_dists(term_decidates>min(years)))-MP(j).Z.termdist(p),'.',...
    %         'color',seas_cmap(zseas(p),:),'markerfacecolor',seas_cmap(zseas(p),:),'linewidth',1,'markersize',10); hold on;
    % end
    set(gca,'xlim',[min(years)-1 max(years)],'ylim',[0,max(term_dists(term_decidates>(min(years)-1)))-min(term_dists(term_decidates>(min(years)-1)))]);
    ylims = get(gca,'ylim'); yticks = get(gca,'ytick');
    set(gca,'ytick',yticks,'yticklabels',yticks/1000,'fontsize',16);
    grid on; drawnow;
    axr.YAxis(1).Color = 'k'; axr.YAxis(2).Color = 'k'; 
    axr.YAxis(2).Label.String = 'Distance from most-retreated terminus (km)'; axr.YAxis(2).Label.FontSize = 16;
    term_leg = legend(pterm,season_names,'Location','northoutside',...
            'Orientation','horizontal');

    %save the figure
    saveas(ts_fig,[root_dir,site_abbrev,'/',site_abbrev,'-terminus-buttressing_',num2str(zcutoff),'m-zcutoff_timeseries.png'],'png'); %save the plots
    clear axr; close(ts_fig);
    % uiwait %advance only after figure is closed


    %add the data to the composite figure
    figure(summaryfig); plot_ind = find(geo_ind == l);
    if ~isempty(plot_ind)
        %navigate to the subplot
        eval(['subplot(subB',num2str(l),');']);
        ax1 = gca;

        %plot the buttressing timeseries
        yyaxis left;
        for k = 1:4
            seas_idx = find(zseas==k);
            for p = 1:length(seas_idx)
                bpM = fill([zdate(seas_idx(p))-0.075,zdate(seas_idx(p))+0.075,zdate(seas_idx(p))+0.075,zdate(seas_idx(p))-0.075,zdate(seas_idx(p))-0.075],...
                    [0,0,butt_Meng(seas_idx(p)),butt_Meng(seas_idx(p)),0]./10^6,seas_cmap(k,:),...
                    'edgecolor','none','facealpha',0.5,'linewidth',1.5); hold on;
            end
        end
        for k = 1:4
            seas_idx = find(zseas==k);
            for p = 1:length(seas_idx)
                bpA = fill([zdate(seas_idx(p))-0.075,zdate(seas_idx(p))+0.075,zdate(seas_idx(p))+0.075,zdate(seas_idx(p))-0.075,zdate(seas_idx(p))-0.075],...
                    [0,0,butt_Amundson(seas_idx(p)),butt_Amundson(seas_idx(p)),0]./10^6,'k',...
                    'edgecolor','none','facealpha',0.5,'linewidth',0.5); hold on;
            end
        end
        set(ax1,'xlim',[min(years)-1 max(years)],'xticklabel',[],'fontsize',14); grid on;
        ax1.YAxis(1).Color = 'k'; set(ax1,'box','on'); ax1.LineWidth = 1.5;
        %add labels
        if plot_locs(plot_ind) >= rows*cols - (cols-1) || (plot_locs(plot_ind) >= 2 && plot_locs(plot_ind) <= 3)
            % xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks); clear xticks;
            if plot_locs(plot_ind) == rows*cols - (cols-1)
                ax1.YAxis(1).Label.String = 'Buttressing (10^6 N/m)';
                ax1.YAxis(1).Label.Position = [2008.4 3.5 -1];
                ax1.YAxis(1).Label.FontSize = 14;
            end
        end
        %crop the y-axis as needed & label bars that extend off the limit
        ylims = get(gca,'ylim'); if max(ylims) > 5; set(gca,'ylim',[0 5]); end
        bM_outlier = find(butt_Meng > 5e6); bA_outlier = find(butt_Amundson > 5e6);
        if ~isempty(bM_outlier)
            if isempty(bA_outlier) %strain-based buttressing data fit in y-axis
                text(zdate(bM_outlier),5.125,num2str(round(butt_Meng(bM_outlier)/10^6,1)),'fontsize',14,'color',seas_cmap(zseas(bM_outlier),:));
                text(zdate(bM_outlier)+1,5.625,'x10^6','fontsize',14,'color','k');
            else %both buttressing estimates exceed y-limits and need to be plotted
                text(zdate(bM_outlier)-1.2,5.625,'(     ,     )x10^6','fontsize',14,'color','k');
                text(zdate(bM_outlier)-0.9,5.5,num2str(round(butt_Meng(bM_outlier)/10^6,0)),'fontsize',14,'color',seas_cmap(zseas(bM_outlier),:));
                text(zdate(bA_outlier)+0.5,5.5,num2str(round(butt_Amundson(bA_outlier)/10^6,0)),'fontsize',14,'color',seas_cmap(zseas(bA_outlier),:));
            end
        end
        clear b*_outlier ylims;

        %plot the terminus timeseries
        yyaxis right; axr = gca;
        plot(term_decidates_sorted,term_pos_sorted,'-k','linewidth',1,'markersize',10); hold on;
        for p = 1:length(term_pos_sorted)
            %add dummy plot for the legend
            if p ==1
                for k = 1:4
                    pterm(k) = plot(term_decidates_sorted(p),NaN,symbol_shape,...
                        'markeredgecolor',seas_cmap(k,:),'markerfacecolor',seas_cmap(k,:),'linewidth',1,'markersize',3); hold on;
                end
            end

            %plot the data
            plot(term_decidates_sorted(p),term_pos_sorted(p),symbol_shape,...
                'color',seas_cmap(term_seas_sorted(p),:),'markerfacecolor',seas_cmap(term_seas_sorted(p),:),'linewidth',1,'markersize',2); hold on;
        end
        % for p = 1:length(MP(j).Z.date)
        %     plot(zdate(p),max(term_dists(term_decidates>min(years)))-MP(j).Z.termdist(p),'.',...
        %         'color',seas_cmap(zseas(p),:),'markerfacecolor',seas_cmap(zseas(p),:),'linewidth',1,'markersize',3); hold on;
        % end
        set(axr,'xlim',[min(years)-1 max(years)],'xticklabel',[],...
            'ylim',[0,max(term_dists(term_decidates>(min(years)-1)))-min(term_dists(term_decidates>(min(years)-1)))]);
        ylims = get(gca,'ylim'); yticks = get(gca,'ytick');
        set(gca,'ytick',yticks,'yticklabels',yticks/1000,'fontsize',14);
        grid on; drawnow;
        axr.YAxis(1).Color = 'k'; axr.YAxis(2).Color = 'k';
        if plot_locs(plot_ind) >= rows*cols - (cols-1)
            xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks); clear xticks;
            set(gca,'yticklabel',yticks/1000);
            xlabel('Year','fontsize',14); 
            if plot_locs(plot_ind) == rows*cols
                axr.YAxis(2).Label.String = 'Distance from most-retreated terminus (km)'; 
                axr.YAxis(2).Label.Position = [2024.4 15e3 -1];
                axr.YAxis(2).Label.FontSize = 14;
            end
        end
        pos = get(gca,'position');
        %add the legend
        if plot_locs(plot_ind) == 2
            term_leg = legend(pterm,season_names,'Location','southoutside',...
                'Orientation','horizontal');
        end
        %shift plot locations
        if mod(plot_locs(plot_ind),3) == 0
            set(gca,'position',[pos(1) pos(2) 0.24 1.1*pos(4)]);
        elseif mod(plot_locs(plot_ind),3) == 2
            set(gca,'position',[pos(1)-0.02 pos(2) 0.24 1.1*pos(4)]);
        else
            set(gca,'position',[pos(1)-0.06 pos(2) 0.24 1.1*pos(4)]);
        end

    end
    drawnow;

    %create plots of biannual terminus position change relative to buttressing
    figure(summaryfig); subdT = subplot(rows,cols,[8,11,14]);
    plot(Bsp(end-length(years)+1:end),dTw(end-length(years)+1:end),symbol_shape,...
        'markeredgecolor','k','linewidth',1,'markerfacecolor',seas_cmap(1,:),'markersize',8); hold on;
    plot(Bsu(end-length(years)+1:end),dTs(end-length(years)+1:end),symbol_shape,...
        'markeredgecolor','k','linewidth',1,'markerfacecolor',seas_cmap(3,:),'markersize',8); hold on;

    clear term*sorted ic;
    clear dts inland_* seaward_idx mean_prof Tdate* term_* Tmos vdist Tyrs zdate* zmos zyrs;
    clear Tdist* zseas Tseas butt_* bp* Havg packing *_profs dVdx dvdist_temp;

end

%label the center plot
figure(summaryfig); subplot(subdT);
%replot the ASG data so it is visible
plot(Bsp(1:length(years)),dTw(1:length(years)),'s',...
    'markeredgecolor','k','linewidth',1,'markerfacecolor',seas_cmap(1,:),'markersize',8); hold on;
plot(Bsu(1:length(years)),dTs(1:length(years)),'s',...
    'markeredgecolor','k','linewidth',1,'markerfacecolor',seas_cmap(3,:),'markersize',8); hold on;
set(subdT,'fontsize',14); grid on; box on;
xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
set(gca,'xlim',[0 2.5e6],'xticklabel',xticks/10^6,'yticklabel',yticks/10^3);
xlabel('Buttressing (10^6 N/m)','fontsize',14); ylabel('Terminus change (km)','fontsize',14); 
% fit only to the most-extensive timeseries (ASG and ZIM)
Bsp_alt = [Bsp(1:length(years)), Bsp(end-length(years)+1:end)]; 
Bsu_alt = [Bsu(1:length(years)), Bsu(end-length(years)+1:end)]; 
dTw_alt = [dTw(1:length(years)), dTw(end-length(years)+1:end)]; 
dTs_alt = [dTs(1:length(years)), dTs(end-length(years)+1:end)];
[fw_alt,gofw_alt] = fit(Bsp_alt(~isnan(Bsp_alt) & ~isnan(dTw_alt)),dTw_alt(~isnan(Bsp_alt) & ~isnan(dTw_alt)),'poly1'); 
[fs_alt,gofs_alt] = fit(Bsu_alt(~isnan(Bsu_alt) & ~isnan(dTs_alt)),dTs_alt(~isnan(Bsu_alt) & ~isnan(dTs_alt)),'poly1'); 
%add extra symbols to call-out ZIM (ommitted from some fits) and ASG (other detailed record)
plot(Bsp(1:length(years)),dTw(1:length(years)),'+w','markersize',6,'linewidth',1); hold on;
plot(Bsu(1:length(years)),dTs(1:length(years)),'+w','markersize',6,'linewidth',1); hold on;
plot(Bsp(end-length(years)+1:end),dTw(end-length(years)+1:end),'xw','markersize',6,'linewidth',1); hold on;
plot(Bsu(end-length(years)+1:end),dTs(end-length(years)+1:end),'xw','markersize',6,'linewidth',1); hold on;
% fit all
[fw,gofw] = fit(Bsp(~isnan(Bsp) & ~isnan(dTw)),dTw(~isnan(Bsp) & ~isnan(dTw)),'poly1'); 
[fs,gofs] = fit(Bsu(~isnan(Bsu) & ~isnan(dTs)),dTs(~isnan(Bsu) & ~isnan(dTs)),'poly1'); 
% % plot the fits
% Bx = linspace(min(Bsp),max(Bsp),10);
% plot(Bx,fw(Bx),'-','color',seas_cmap(1,:),'linewidth',2); plot(Bx,fs(Bx),'-','color',seas_cmap(3,:),'linewidth',2); hold on;
% plot(Bx,fw_alt(Bx),'--','color',seas_cmap(1,:),'linewidth',2); plot(Bx,fs_alt(Bx),'--','color',seas_cmap(3,:),'linewidth',2); hold on;
set(gca,'ytick',[-8*10^3:2*10^3:8*10^3],'yticklabel',[-8:2:8]);
ylims = get(gca,'ylim'); xlims = get(gca,'xlim');
%adjust subplot position
subplot(subB11); pos = get(gca,'position');
subplot(subdT); subpos = get(gca,'position'); set(gca,'linewidth',1.5);
set(subdT,'position',[pos(1) subpos(2)-0.06 pos(3) subpos(4)+0.12]);
% add legend
rectangle('position',[min(xlims)+0.2*range(xlims), min(ylims)+0.875*range(ylims),0.775*range(xlims),0.1*range(ylims)],'facecolor', ...
    'w','edgecolor','k','linewidth',1.5);
%symbols: big3 vs others
plot(min(xlims)+0.23*range(xlims), min(ylims)+0.945*range(ylims),...
    'd','markeredgecolor','k','markerfacecolor','w','linewidth',1,'markersize',8);
plot(min(xlims)+0.23*range(xlims), min(ylims)+0.905*range(ylims),...
    's','markeredgecolor','k','markerfacecolor','w','linewidth',1,'markersize',8);
text(min(xlims)+0.265*range(xlims),min(ylims)+0.95*range(ylims),['D>15 Gt yr^{-1}'],'fontsize',14);
text(min(xlims)+0.265*range(xlims),min(ylims)+0.91*range(ylims),['D<15 Gt yr^{-1}'],'fontsize',14);
%colors: winter vs summer
rectangle('position',[min(xlims)+0.65*range(xlims), min(ylims)+0.935*range(ylims),0.04*range(xlims),0.02*range(ylims)],'facecolor', ...
    seas_cmap(1,:),'edgecolor','none','linewidth',1);
rectangle('position',[min(xlims)+0.65*range(xlims), min(ylims)+0.895*range(ylims),0.04*range(xlims),0.02*range(ylims)],'facecolor', ...
    seas_cmap(3,:),'edgecolor','none','linewidth',1);
% plot(min(xlims)+0.725*range(xlims), min(ylims)+0.955*range(ylims),...
%     symbol_shape,'markeredgecolor','k','markerfacecolor',seas_cmap(1,:),'linewidth',1);
% plot(min(xlims)+0.725*range(xlims), min(ylims)+0.925*range(ylims),...
%     symbol_shape,'markeredgecolor','k','markerfacecolor',seas_cmap(3,:),'linewidth',1);
text(min(xlims)+0.71*range(xlims),min(ylims)+0.95*range(ylims),['winter'],'fontsize',14);
text(min(xlims)+0.71*range(xlims),min(ylims)+0.91*range(ylims),['summer'],'fontsize',14);
%ASG & ZIM
rectangle('position',[min(xlims)+0.38*range(xlims), min(ylims)+0.015*range(ylims),0.60*range(xlims),0.09*range(ylims)],'facecolor', ...
    'w','edgecolor','k','linewidth',1.5);
plot(min(xlims)+0.415*range(xlims), min(ylims)+0.08*range(ylims),...
    's','markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5],'linewidth',1,'markersize',8);
plot(min(xlims)+0.415*range(xlims), min(ylims)+0.08*range(ylims),'+w','markersize',6,'linewidth',1);
plot(min(xlims)+0.415*range(xlims), min(ylims)+0.04*range(ylims),...
    'd','markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5],'linewidth',1,'markersize',8);
plot(min(xlims)+0.415*range(xlims), min(ylims)+0.04*range(ylims),'xw','markersize',6,'linewidth',1);
text(min(xlims)+0.46*range(xlims),min(ylims)+0.08*range(ylims),['Nunatakassaap'],'fontsize',14);
text(min(xlims)+0.46*range(xlims),min(ylims)+0.04*range(ylims),['Zachari',char(230),' Isstr',char(248),'m'],'fontsize',14);
% % add legend for fits
% rectangle('position',[min(xlims)+0.2*range(xlims), min(ylims)+0.025*range(ylims),0.775*range(xlims),0.13*range(ylims)],'facecolor', ...
%     'w','edgecolor','k','linewidth',1.5);
% %all data
% plot([min(xlims)+0.225*range(xlims),min(xlims)+0.3*range(xlims)], [min(ylims)+0.12*range(ylims),min(ylims)+0.12*range(ylims)],...
%     '-','color',seas_cmap(1,:),'linewidth',2);
% plot([min(xlims)+0.225*range(xlims),min(xlims)+0.3*range(xlims)], [min(ylims)+0.055*range(ylims),min(ylims)+0.055*range(ylims)],...
%     '-','color',seas_cmap(3,:),'linewidth',2);
% text(min(xlims)+0.325*range(xlims),min(ylims)+0.12*range(ylims),['R^2_{all}=',num2str(round(gofw.rsquare,2))]);
% text(min(xlims)+0.325*range(xlims),min(ylims)+0.055*range(ylims),['R^2_{all}=',num2str(round(gofs.rsquare,2))]);
% %ASG & ZIM only
% plot([min(xlims)+0.575*range(xlims),min(xlims)+0.675*range(xlims)], [min(ylims)+0.12*range(ylims),min(ylims)+0.12*range(ylims)],...
%     '--','color',seas_cmap(1,:),'linewidth',2);
% plot([min(xlims)+0.575*range(xlims),min(xlims)+0.675*range(xlims)], [min(ylims)+0.055*range(ylims),min(ylims)+0.055*range(ylims)],...
%     '--','color',seas_cmap(3,:),'linewidth',2);
% text(min(xlims)+0.70*range(xlims),min(ylims)+0.12*range(ylims),['R^2_{decadal}=',num2str(round(gofw_alt.rsquare,2))]);
% text(min(xlims)+0.70*range(xlims),min(ylims)+0.055*range(ylims),['R^2_{decadal}=',num2str(round(gofs_alt.rsquare,2))]);
disp(['Winter advance = ',num2str(fw.p1*10^6/1000), 'km for each 10^6 N/m increase in buttressing']);
disp(['Winter fit RMSE = ',num2str(gofw.rmse/1000), ' km and R^2 = ',num2str(gofw.rsquare)]);
disp(['Summer retreat = ',num2str(fs.p1*10^6/1000), 'km for each 10^6 N/m increase in buttressing']);
disp(['Summer fit RMSE = ',num2str(gofs.rmse/1000), ' km and R^2 = ',num2str(gofs.rsquare)]);

%save the GrIS-wide timeseries
saveas(summaryfig,[root_dir,'GrIS-melange_terminus-buttressing_',num2str(zcutoff),'m-zcutoff_timeseries.png'],'png');
exportgraphics(summaryfig,[root_dir,'GrIS-melange_terminus-buttressing_',num2str(zcutoff),'m-zcutoff_timeseries.tif'],Resolution=600);

%add labels to the figure comparing buttressing and terminus position
figure(termfig);
subplot(sub1); set(sub1,'fontsize',16);
xlabel('Thickness (m)','fontsize',16); ylabel('Buttressing (N/m)','fontsize',16); 
subplot(sub2); set(sub2,'fontsize',16);
xlabel('Buttressing (N/m)','fontsize',16); ylabel('Terminus change rate (m/d)','fontsize',16); 
% close(termfig);


% create violin plots to show seasonal profile distributions
prof_fig = figure; set(prof_fig,'position',[350 50 900 1200]);
Hdist_cat = categorical(Hprofs_dist'); max_dist = 12000;
% convert thickness data to vectors for plotting
dist_groups = []; H_vals = []; season_groups = []; dataflag = [];
for k = 2:4
    %isolate the data for the season
    ygroupdata = profs(k).H';
    for l = 1:size(ygroupdata,1)
        obsno(l) = sum(~isnan(ygroupdata(l,:)));
    end
    % endref = find(obsno>1,1,'last');
    endref = find(Hprofs_dist<=max_dist,1,'last');
    if ~isempty(endref)
        dataflag = [dataflag; k];
        ygroupdata = ygroupdata(1:endref,:);
        xgroupdata = categorical(repmat(Hprofs_dist(1:endref)'/1000,1,size(ygroupdata,2)));
        
        %convert to column vectors
        xgroupcolumn = xgroupdata(~isnan(ygroupdata)); ycolumn = ygroupdata(~isnan(ygroupdata));
        cgroupdata = categorical(repelem(season_names(k),length(xgroupcolumn)));
        
        %concatenate
        dist_groups = [dist_groups; xgroupcolumn];
        H_vals = [H_vals; ycolumn];
        season_groups = [season_groups; cgroupdata'];
    end

    %plot a bar plot of the observation count profile
    subplot(3,3,k-1);
    bp = bar(Hprofs_dist(1:endref)'/1000,obsno(1:endref));
    set(gca,'ylim',[0,150],'xlim',[-1,max_dist/1000+1],'xticklabel',[],'fontsize',16); 
    if k == 2; ylabel('Count','fontsize',16); else; set(gca,'yticklabel',[]); end
    bp.FaceColor = seas_cmap(k,:); bp.FaceAlpha = 0.5; clear bp;
    pos = get(gca,'position'); set(gca,'position',[pos(1)-(4-k)*0.03 0.86 1.35*pos(3) 0.1]);
    ylims = get(gca,'ylim'); text(max_dist/1000,0.8*max(ylims),[alphabet(k-1),') '],'fontsize',16);

    %plot individual season data in its own panel
    subplot(3,3,3+k-1);
    vp = violinplot(xgroupcolumn,ycolumn); hold on;
    vp.FaceColor = seas_cmap(k,:); vp.FaceAlpha = 0.5; clear vp;
    plot(categorical(Hprofs_dist(1:endref)'/1000),nanmedian(ygroupdata,2),'*k');
    set(gca,'ylim',[0,150],'fontsize',16,'xticklabel',[]); box on; ylims = get(gca,'ylim');
    pos = get(gca,'position'); set(gca,'position',[pos(1)-(4-k)*0.03 0.48 1.35*pos(3) 0.35]);
    if k == 2
        ylabel('Mean thickness (m)','fontsize',16); 
    else
        set(gca,'yticklabel',[]);
    end
    text(max_dist/1000+1,0.95*max(ylims),[alphabet(3+k-1),') '],'fontsize',16);
    % %add number of observations
    % for l = 1:endref
    %     if obsno(l) > 100
    %         text(l,0.92*max(ylims),['(',num2str(obsno(l)),')'],'fontsize',12,'Rotation',90);
    %     else
    %         text(l,0.92*max(ylims),['(',num2str(obsno(l)),')'],'fontsize',12,'Rotation',90);
    %     end
    % end

    clear ygroupdata xgroup* ycolumn cgroupdata; clear obsno;
end
% plot all the data in a single panel
% vp = violinplot(dist_groups,H_vals,GroupByColor=season_groups);
% for k = 1:length(vp)
%     vp(k).FaceColor = seas_cmap(dataflag(k),:);
% end
% set(gca,'ylim',[0,150],'fontsize',16); box on;

% convert packing density data to vectors for plotting
dist_groups = []; p_vals = []; season_groups = []; dataflag = [];
for k = 2:4
    %isolate the data for the season
    ygroupdata = 100*profs(k).p';
    for l = 1:size(ygroupdata,1)
        obsno(l) = sum(~isnan(ygroupdata(l,:)));
    end
    % endref = find(obsno>1,1,'last');
    endref = find(Hprofs_dist<=max_dist,1,'last');
    if ~isempty(endref)
        dataflag = [dataflag; k];
        ygroupdata = ygroupdata(1:endref,:);
        xgroupdata = categorical(repmat(Hprofs_dist(1:endref)'/1000,1,size(ygroupdata,2)));
        
        %convert to column vectors
        xgroupcolumn = xgroupdata(~isnan(ygroupdata)); ycolumn = ygroupdata(~isnan(ygroupdata));
        cgroupdata = categorical(repelem(season_names(k),length(xgroupcolumn)));

        %concatenate
        dist_groups = [dist_groups; xgroupcolumn];
        p_vals = [p_vals; ycolumn];
        season_groups = [season_groups; cgroupdata'];
    end
    clear obsno;

    %plot individual season data in its own panel
    subplot(3,3,6+k-1);
    vp = violinplot(xgroupcolumn,ycolumn); hold on;
    vp.FaceColor = seas_cmap(k,:); vp.FaceAlpha = 0.5; clear vp;
    plot(categorical(Hprofs_dist(1:endref)'/1000),nanmedian(ygroupdata,2),'*k');
    set(gca,'ylim',[0,100],'fontsize',16); box on; ylims = get(gca,'ylim');
    pos = get(gca,'position'); set(gca,'position',[pos(1)-(4-k)*0.03 0.10 1.35*pos(3) 0.35]);
    xlabel('Distance from terminus (km)','fontsize',16);
    if k == 2
        ylabel('Packing density (%)','fontsize',16);
    else
        set(gca,'yticklabel',[]);
    end
    text(max_dist/1000+1,0.95*max(ylims),[alphabet(6+k-1),') '],'fontsize',16);

    clear ygroupdata xgroup* ycolumn cgroupdata;
end

% % convert strain rate data to vectors for plotting
% dist_groups = []; V_vals = []; season_groups = []; dataflag = [];
% for k = 2:4
%     %isolate the data for the season
%     ygroupdata = profs(k).V';
% 
%     %truncate the along-fjord plot extent
%     for l = 1:size(ygroupdata,1)
%         obsno(l) = sum(~isnan(ygroupdata(l,:)));
%     end
%     % endref = find(obsno>1,1,'last');
%     endref = find(Hprofs_dist<=max_dist,1,'last');
% 
%     %identify the range of values & fit kernel density functions within
%     % %that range to the data at that location
%     % for l = 1:size(ygroupdata,1)
%     %     pts = linspace(min(ygroupdata(l,:),[],'omitmissing'), max(ygroupdata(l,:),[],'omitmissing'), size(ygroupdata,2));
%     %     [densVals, evalPts] = ksdensity(ygroupdata(l,:), pts);
%     % end
% 
%     if ~isempty(endref)
%         dataflag = [dataflag; k];
%         ygroupdata = ygroupdata(1:endref,:);
%         xgroupdata = categorical(repmat(Hprofs_dist(1:endref)'/1000,1,size(ygroupdata,2)));
% 
%         %convert to column vectors
%         xgroupcolumn = xgroupdata(~isnan(ygroupdata)); ycolumn = ygroupdata(~isnan(ygroupdata));
%         cgroupdata = categorical(repelem(season_names(k),length(xgroupcolumn)));
% 
%         %concatenate
%         dist_groups = [dist_groups; xgroupcolumn];
%         V_vals = [V_vals; ycolumn];
%         season_groups = [season_groups; cgroupdata'];
%     end
%     clear obsno;
% 
%     %plot individual season data in its own panel
%     subplot(3,3,6+k-1);
%     vp = violinplot(xgroupcolumn,ycolumn);
%     vp.FaceColor = seas_cmap(k,:); vp.FaceAlpha = 0.5; clear vp;
%     set(gca,'ylim',[0,2e4],'fontsize',16); box on;
%     pos = get(gca,'position'); set(gca,'position',[pos(1)-(4-k)*0.03 pos(2)-0.05 1.35*pos(3) 1.2*pos(4)]);
% 
%     clear ygroupdata xgroup* ycolumn cgroupdata;
% end

%save the violin plot profiles
saveas(prof_fig,[root_dir,'GrIS-melange_',num2str(zcutoff),'m-zcutoff_violin-profiles.png'],'png');
exportgraphics(prof_fig,[root_dir,'GrIS-melange_',num2str(zcutoff),'m-zcutoff_violin-profiles.tif'],Resolution=600);


% close all;



%% create terminus position plots
load([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat']);

for l = 1:length(sitenames)
    %identify the index in the structure for the specified sitename
    j = find(strcmp(string(saved_sites),sitenames(l,:)));

    %specify site plotting conventions
    site_abbrev = MP(j).name; disp(site_abbrev);
    close all; drawnow;
    ts_fig = figure; set(ts_fig,'position',[50 50 600 1200]);

    %extract the transect_inc from the files (only works if only one
    %version is saved for each site but it's needed to account for adaptive
    %transect spacing depending on site size)
    shp_files = dir([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'*.shp']);
    for k = 1:length(shp_files)
        if contains(shp_files(k).name,['transects_'])
            transect_inc = str2num(shp_files(k).name(end-8:end-5));
        end
    end

    %load the centerline & create a distance vector
    C = shaperead([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'-centerline.shp']);
    % centerline_dist = [0, cumsum(sqrt((C.X(2:end)-C.X(1:end-1)).^2 + (C.Y(2:end)-C.Y(1:end-1)).^2))]';
    % for l = 1:length(MP(j).V.X)
    %     dists = sqrt((MP(j).V.X(l)-C.X).^2 + (MP(j).V.Y(l)-C.Y).^2);
    %     tran_dist(l) = centerline_dist(find(dists == min(dists)));
    %     clear dists;
    % end

    %DEM data
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        zdatetime(p) = datetime(MP(j).Z.date{p},'InputFormat','yyyyMMdd');
        zyrs(p) = year(zdatetime(p)); zmos(p) = month(zdatetime(p));
        for k = 1:4
            if ismember(zmos(p),seasons(k,:))==1
                zseas(p) = k;
            end
        end
    end
    mean_prof = nanmean(MP(j).Z.transectZavg,2);
    seaward_idx = find(~isnan(mean_prof)==1,1,'first');
    inland_idx = find(~isnan(mean_prof)==1,1,'last');

    %intersect the centerline with each time-stamped melange outline
    for p = 1:length(MP(j).T.date)
        try 
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        catch
            MP(j).T.date{p} = string(str2num(MP(j).T.date{p})-1);
            Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        end
        Tdate(p) = convert_to_decimaldate(char(MP(j).T.date{p}));
        Tdatetime(p) = datetime(MP(j).T.date{p},'InputFormat','yyyyMMdd');
        Tyrs(p) = year(Tdatetime(p)); Tmos(p) = month(Tdatetime(p));
        dists = sqrt((MP(j).T.termX(1,p)-C.X).^2 + (MP(j).T.termY(1,p)-C.Y).^2);
        [~,term_ind] = min(dists);
        % MP(j).T.termdist(1,p) = centerline_dist(term_ind);
        for k = 1:4
            if ismember(Tmos(p),seasons(k,:))==1
                Tseas(p) = k;
            end
        end
        clear dists term_ind;
    end
    MP(j).T.termdist(Tyrs < min(years) | Tyrs > max(years)) = NaN;

    %compile the data for identifying terminus gaps
    term_decidates = [Tdate, zdate(term_trace ==1)]; term_dates = [Tdatetime, zdatetime(term_trace ==1)]; 
    term_X = [MP(j).T.termX, MP(j).Z.termX(term_trace ==1)];
    term_Y = [MP(j).T.termY, MP(j).Z.termY(term_trace ==1)];
    term_dists = [MP(j).T.termdist, MP(j).Z.termdist(term_trace ==1)];

    %plot the image-based termini first, then DEM
    for p = 1:length(MP(j).T.date)
        %add dummy plot for the legend
        if p ==1
            for k = 1:4
            pterm(k) = plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).T.termdist(p),Tdate(p),'x',...
            'color',seas_cmap(k,:),'linewidth',2,'markersize',10); hold on;
            end
        end

        %plot the data
        plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).T.termdist(p),Tdate(p),'x',...
            'color',seas_cmap(Tseas(p),:),'linewidth',2,'markersize',10); hold on;
    end
    for p = 1:length(MP(j).Z.date)
        if term_trace(p) == 0 && isnan(MP(j).Z.termdist(p))
            disp(['Still need terminus data for ',MP(j).name,' ',char(MP(j).Z.date(p))]);
        else
            plot(max([MP(j).T.termdist,MP(j).Z.termdist])-MP(j).Z.termdist(p),zdate(p),'s',...
                'color',seas_cmap(zseas(p),:),'linewidth',1,'markerfacecolor','none','markersize',10); hold on;
        end
    end
    set(gca,'xlim',[0,max([MP(j).T.termdist,MP(j).Z.termdist])-min([MP(j).T.termdist,MP(j).Z.termdist])],'ylim',[min(years) max(years)]);
    xlims = get(gca,'xlim'); xticks = get(gca,'xtick');
    set(gca,'xtick',xticks,'xticklabels',xticks/1000,'fontsize',16);
    term_leg = legend(pterm,season_names,'Location','northoutside',...
            'Orientation','horizontal');
    xlabel('Distance from most-retreated terminus (km)','fontsize',16); ylabel('Year','fontsize',16);
    grid on; drawnow;
    
    %save the figure & updated terminus positions
    % save([root_dir,'GrIS-melange-characteristics.mat'],'MP','-v7.3');
    saveas(ts_fig,[root_dir,MP(j).name,'/',MP(j).name,'-terminus-timeseries-plot.png'],'png'); %save the plots
    % uiwait %advance only after figure is closed

    %create characteristic seasonal curves for terminus position
    MP(j).T.termanom = NaN(size(MP(j).T.termdist));
    for p = 1:length(years)
        yr_idx = find(Tyrs == years(p));
        if ~isempty(yr_idx)
            %identify the months in the specificyear
            mos_yr = Tmos(yr_idx);

            %calculate seasonal medians
            for k = 1:4
                Tdist_seas(k,p) = nanmedian(MP(j).T.termdist(yr_idx(ismember(mos_yr,seasons(k,:))==1)));
            end

            %calculate the mean of the seasonal positions
            Tmean = nanmean(Tdist_seas(:,p));

            %subtract seasonal medians from the annual median so that
            %positive anomalies indicate terminus advance (farther from
            %most retreated terminus) & negative anomalies indicate retreat
            MP(j).T.termanom(yr_idx) = Tmean - MP(j).T.termdist(yr_idx);
            Tdist_seasanom(:,p) = Tmean - Tdist_seas(:,p);
            
            clear Tmean mos_y;
        end
        clear yr_idx;
    end

    %calculate monthly anomalies
    MP(j).T.termanom_mo = NaN(12,1); MP(j).T.termanom_morange = NaN(12,2); 
    for k = 1:12
        MP(j).T.termanom_mo(k,1) = nanmean(MP(j).T.termanom(ismember(Tmos,k)==1));
        MP(j).T.termanom_morange(k,:) = prctile(MP(j).T.termanom(ismember(Tmos,k)==1),[5,95],2);
    end

    clear C centerline_dist dts inland_idx seaward_idx mean_prof Tdate* term_* Tmos tran_dist Tyrs zdate* zmos zyrs;
    clear Tdist* zseas Tseas;
end
%resave the data (as needed)
save([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat'],'MP','-v7.3');



%% plot all the characteristic monthly seasonal terminus anomaly timeseries
seasterm_fig = figure; set(seasterm_fig,'position',[50 50 1200 400]);
sub1 = subplot(3,1,1); sub2 = subplot(3,1,2); sub3 = subplot(3,1,3);
% term_cmap = [223,194,125; 125,125,125]./255; %shades of brown & gray
% term_cmap = colormap(gray(length(MP)));
term_cmap = cmocean('phase',15);

%normalize the data or keep as full anomalies?
normterm = 0; %1 = normalize, 0 = raw

%fill the background with seasonal colors
for k = 1:3
    eval(['subplot(sub',num2str(k),');'])
    if normterm == 1; yscale = 1; else; yscale = 2000; end
    fill([seasons(1,1),13,13,seasons(1,1),seasons(1,1)],yscale*[-1,-1,1,1,-1],seas_cmap(1,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    fill([0,seasons(2,1),seasons(2,1),0,0],yscale*[-1,-1,1,1,-1],seas_cmap(1,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    fill([seasons(2,1),seasons(3,1),seasons(3,1),seasons(2,1),seasons(2,1)],yscale*[-1,-1,1,1,-1],seas_cmap(2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    fill([seasons(3,1),seasons(4,1),seasons(4,1),seasons(3,1),seasons(3,1)],yscale*[-1,-1,1,1,-1],seas_cmap(3,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    fill([seasons(4,1),seasons(1,1),seasons(1,1),seasons(4,1),seasons(4,1)],yscale*[-1,-1,1,1,-1],seas_cmap(4,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    set(gca,'xlim',[1,13],'xtick',[],'fontsize',16); grid on;
    if k == 2
        ylabel('Normalized seasonal terminus anomaly','fontsize',16);
    elseif k == 3
        xlabel('Month','fontsize',16); 
        set(gca,'xtick',[1:12]);
    end
end
%plot the data on top
pl_ref = 1;
for j = 1:length(geo_ind)
    if ~contains(MP(geo_ind(j)).name,'KBG')
        disp(MP(geo_ind(j)).name)
        %decide the color for the line based on when it starts to retreat
        if nanmean(MP(geo_ind(j)).T.termanom_mo(4:5)) > 30+MP(geo_ind(j)).T.termanom_mo(6)
            cmap_ind = 1; %early retreater!
        else
            if MP(geo_ind(j)).T.termanom_mo(6) > 30+MP(geo_ind(j)).T.termanom_mo(7)
                cmap_ind = 2; %June peak
            else
                cmap_ind = 3; %July peak
            end
        end
        dbutt = round(100*nanmean(squeeze(MP(geo_ind(j)).B.butt_Meng(1,3,:)))./nanmean(squeeze(MP(geo_ind(j)).B.butt_Meng(1,2,:))));
        if isnan(dbutt); dbutt = 1; end
        if dbutt > 100; lstyle = '-'; else; lstyle = '--'; end

        %normalize as needed
        if normterm == 1
            tanoms = MP(geo_ind(j)).T.termanom_mo./max(abs(MP(geo_ind(j)).T.termanom_mo));
        else
            tanoms = MP(geo_ind(j)).T.termanom_mo;
        end


        % plot based on site size & timing of seasonal retreat
        if ismember(MP(geo_ind(j)).name,big3)
            eval(['subplot(sub',num2str(cmap_ind),');'])
            pl(pl_ref) = plot([1:12]+0.5,tanoms,'linestyle',lstyle,'color',term_cmap(j,:),'linewidth',2); hold on;
                plot([1:12]+0.5,tanoms,'d',...
                    'color',term_cmap(j,:),'linewidth',1,'markerfacecolor',term_cmap(j,:)); hold on;
        else
            eval(['subplot(sub',num2str(cmap_ind),');'])
            pl(pl_ref) = plot([1:12]+0.5,tanoms,'linestyle',lstyle,'color',term_cmap(j,:),'linewidth',2); hold on;
                plot([1:12]+0.5,tanoms,'s',...
                    'color',term_cmap(j,:),'linewidth',1,'markerfacecolor',term_cmap(j,:)); hold on;
        end
        % site_ind(pl_ref) = strmatch(site_abbrev,geo_order);  
        site_ind(pl_ref) = geo_ind(j);  pl_ref = pl_ref+1; clear tanoms;
    end
end
term_leg = legend(pl,geo_names); set(term_leg,'location','eastoutside'); 
sub1pos = get(sub1,'position'); sub2pos = get(sub2,'position'); sub3pos = get(sub3,'position'); 
set(sub2,'position',[sub2pos(1) sub2pos(2) 0.63 1.1*sub2pos(4)]);
set(sub3,'position',[sub3pos(1) sub3pos(2) 0.63 1.1*sub3pos(4)]);
%format y-limits
for k = 1:3
    eval(['subplot(sub',num2str(k),');']);
    if normterm == 0
        ylims = get(gca,'ylim');
        if range(ylims) > 3000
            set(gca,'ylim',[1000*floor(min(ylims)/1000),1000*ceil(max(ylims)/1000)],'ytick',[1000*floor(min(ylims)/1000):1000:1000*ceil(max(ylims)/1000)]);
        else
            set(gca,'ylim',[500*floor(min(ylims)/500),500*ceil(max(ylims)/500)],'ytick',[500*floor(min(ylims)/500):500:500*ceil(max(ylims)/500)]);
        end
    end

    %add labels
    if k == 1
        text(1.1,min(get(gca,'ylim'))+0.15*range(get(gca,'ylim')),'a) April/May terminus maximum','fontsize',14); 
    elseif k == 2
        subplot(sub2); text(1.1,min(get(gca,'ylim'))+0.15*range(get(gca,'ylim')),'b) June terminus maximum','fontsize',14); 
    else
        subplot(sub3); text(1.1,min(get(gca,'ylim'))+0.15*range(get(gca,'ylim')),'c) July terminus maximum','fontsize',14);
    end
end
legpos = get(term_leg,'position'); set(term_leg,'position',[legpos(1) 0.2 legpos(3) legpos(4)]);
set(sub1,'position',[sub1pos(1) sub1pos(2) 0.63 1.1*sub1pos(4)]);
if normterm == 1
    saveas(seasterm_fig,[root_dir,'GrIS-terminus-seasonal-normanomalies_plot.png'],'png'); %save the plots
    exportgraphics(seasterm_fig,[root_dir,'GrIS-terminus-seasonal-normanomalies_plot.tif'],Resolution=600);
else
    saveas(seasterm_fig,[root_dir,'GrIS-terminus-seasonal-anomalies_plot.png'],'png'); %save the plots
    exportgraphics(seasterm_fig,[root_dir,'GrIS-terminus-seasonal-anomalies_plot.tif'],Resolution=600);
end
close all; 
clear pl term_cmap site_ind;
disp('Done plotting terminus timeseries');



%% create overview maps for each site
if zcutoff == 3
cd(root_dir);
load([root_dir,'GrIS-melange-characteristics_',num2str(zcutoff),'m-zcutoff.mat']);

%plot option
suffix = 'site-map_bigfont'; %add 'bigfont' to the suffix if using 30pt
map_font = 28; %default font = 16, bigfont = 28
if contains(suffix,'bigfont')
    x_shift = 0.05;
else
    x_shift = 0.02;
end

%iterate
for l = 1:length(sitenames)
    close all; drawnow;

    %identify the index in the structure for the specified sitename
    j = find(strcmp(string(saved_sites),sitenames(l,:)));

    %specify site plotting conventions
    site_abbrev = MP(j).name; disp(site_abbrev);
    output_dir = [root_dir,site_abbrev,'/'];
    
    %load the melange masks
    cd([root_dir,site_abbrev]);
    load([site_abbrev,'-melange-masks.mat']); %load the melange mask file

    %load the reference image
    ims = dir([root_dir,site_abbrev,'/images/',site_abbrev,'*S2*.tif']); 
    ref_image = [ims(1).folder,'/',ims(1).name];
    clear ims;
    %load the image geotiff
    [I,R] = readgeoraster(ref_image);
    im.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
    im.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
    if size(I,3) == 1 %panchromatic
        im.z = double(I);
    elseif size(I,3) == 3 %RGB
        im.z = double(I);
    else %all bands
        im.z = double(I(:,:,2:4));
    end
    %crop the image to adjust brightnesses appropriately
    xinds = [find(im.x<=min(melmask.uncropped.x),1,'last')-500/R.CellExtentInWorldX, find(im.x<=max(melmask.uncropped.x),1,'last')+500/R.CellExtentInWorldX];
    yinds = [find(im.y>=max(melmask.uncropped.y),1,'last')-500/R.CellExtentInWorldY, find(im.y<=min(melmask.uncropped.y),1,'first')+500/R.CellExtentInWorldY];
    xinds(xinds<1) = 1; xinds(xinds>length(im.x)) = length(im.x);
    yinds(yinds<1) = 1; yinds(yinds>length(im.y)) = length(im.y);
    im_subset = im.z(min(yinds):max(yinds),min(xinds):max(xinds),:);
    im_subset = im_subset./max(max(im_subset));
    clear I R;

    %isolate terminus positions from the dated melange masks
    for p = 1:length(melmask.dated)
        [~,on] = inpolygon(melmask.dated(p).x,melmask.dated(p).y,melmask.uncropped.x,melmask.uncropped.y);
        mask_xrange(p,:) = [min(melmask.dated(p).x),max(melmask.dated(p).x)]; 
        mask_yrange(p,:) = [min(melmask.dated(p).y),max(melmask.dated(p).y)];
        melmask.dated(p).x(on) = []; melmask.dated(p).y(on) = [];
    end
    %grab info about the DEMs & corresponding terminus positions
    for p = 1:length(MP(j).Z.date)
        term_trace(p) = MP(j).Z.termflag(p);
        zdate(p) = convert_to_decimaldate(char(MP(j).Z.date(p)));
        datest = char(MP(j).Z.date(p));
        zyrs(p) = str2num(datest(1:4)); zmos(p) = str2num(datest(5:6));
        for l = 1:4
            if ismember(zmos(p),seasons(l,:))
                zseas(p) = l;
            end
        end
        clear datest;
    end

    %add manual delineations from images
    cd([root_dir,site_abbrev,'/termini/']);
    % termfiles = dir('*term*.shp');
    if contains(site_abbrev,'KBG') %skip looking for manual delineations other than termpicks for KBG
        termfile_names = {[site_abbrev,'_termpicks.shp']}; %specify consistent data sources
    else
        termfile_names = {[site_abbrev,'_termpicks.shp'],[site_abbrev,'-termini-EPSG3413_2019-2023.shp']}; %specify consistent data sources
    end
    T_inds = []; YYYYMMDD = []; Tseas = []; Tyrs = []; Tmos = []; TX = []; TY = [];
    for m = 1:length(termfile_names)
        % load the shapefile
        % if ismember(termfiles(l).name,termfile_names)
        term = shaperead(char(termfile_names(m))); proj_check = 0;
        for k = 1:length(term)
            term_date(k) = datenum(term(k).Date,'yyyy-mm-dd');
            [xis,~] = polyxpoly(term(k).X,term(k).Y,MP(j).V.X,MP(j).V.Y);
            if isempty(xis); proj_check = proj_check+1; end
        end
        %reproject data as needed (assuming conversion from EPSG4326 > EPSG3413)
        if proj_check == length(term)
            for k = 1:length(term)
                [x,y] = wgs2ps(term(k).X,term(k).Y) ;
                term(k).X = []; term(k).Y = [];
                term(k).X = x; term(k).Y = y; clear x y;
            end
        end
        [~,idx] = sort(term_date);
        for k = 1:length(term)
            sorted_term(k) = term(idx(k));
        end
        clear term; term = sorted_term; clear sorted_term term_date;
        %convert format of TermPicks dates
        bad_ind = [];
        for k = 1:length(term)
            %isolate the year, moth, and day from the date
            if ~isfield(term(k),'Year')
                term(k).Year = num2str(year(term(k).Date));
                term(k).Month = num2str(month(term(k).Date));
                term(k).Day = num2str(day(term(k).Date));
            else
                if k == 1 %wipe out existing data in case format is inconsistent with what I want
                    term = rmfield(term,{'Year','Month','Day'});
                end
                term(k).Year = num2str(year(term(k).Date));
                term(k).Month = num2str(month(term(k).Date));
                term(k).Day = num2str(day(term(k).Date));
            end

            %convert date to same format as used for elevation data
            if length(term(k).Month) == 1; term(k).Month = ['0',term(k).Month]; end
            if length(term(k).Day) == 1; term(k).Day = ['0',term(k).Day]; end
            YYYYMMDD(k) = string([term(k).Year,term(k).Month,term(k).Day]);
            %convert to datetimes
            try
                datestr(k) = datetime(term(k).Date,'InputFormat','yyyy-mm-dd');
            catch
                bad_ind = [bad_ind,k];
            end
        end
        %remove data with erroneous dates
        term(bad_ind) = [];
        YYYYMMDD(bad_ind) = [];

        %concatenate data from different datasets
        if m == 1
            cat_ref = 0;
        else
            cat_ref = length(T);
        end
        for k = 1:length(term)
            T(cat_ref+k).X = term(k).X; T(cat_ref+k).Y = term(k).Y;
            Tyrs(cat_ref+k,1) = str2num(term(k).Year);
            Tmos(cat_ref+k,1) = str2num(term(k).Month);
            for seas = 1:4
                if ismember(str2num(term(k).Month),seasons(seas,:))
                    Tseas(cat_ref+k,1) = seas;
                end
            end
        end
        %clear variables
        clear term YYYYMMDD datestr idx seas;
        % end
    end


    %identify the centerline points from which data were extracted and
    %create a colormap for the bins
    Dsubs = dir([root_dir,site_abbrev,'/',site_abbrev,'*-iceberg-distribution-subsets.csv']);
    size_classes = [];
    for p = 1:length(Dsubs)
        D = readtable([root_dir,site_abbrev,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        name_split = split(Dsubs(p).name,'-',2);
        berg_dates(p) = datetime(name_split{2},'Inputformat','yyyyMMdd'); clear name_split;

        %identify observational limits along the centerline
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes(p,:) = sum(~isnan(berg_nos),1);
        seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
        inland_ext(p) = find(size_classes(p,:)>0,1,'last');
        clear D;
    end
    clear Dsubs;

    %create a colormap for illustrating melange bins
    % tran_cmap = cmocean('deep',(max(inland_ext)-min(seaward_ext)+1));
    yellow_end = [255,255,178]/255; red_end = [189,0,38]/255; colorgrad = (yellow_end-red_end)/(max(inland_ext)-min(seaward_ext)+1);
    for k = 1:(max(inland_ext)-min(seaward_ext)+1)
        tran_cmap(k,:) = yellow_end-(k-1)*colorgrad;
    end
    tran_cmap = flipud(tran_cmap);

    %create an overview map
    map_fig = figure; set(map_fig,'position',[850 50 900 800]); ax1 = axes;
    % imagesc(ax1,im.x(min(xinds):max(xinds)),im.y(min(yinds):max(yinds)),imadjust(im_subset,[],[])); axis xy equal; 
    imagesc(ax1,im.x(min(xinds):max(xinds)),im.y(min(yinds):max(yinds)),im_subset); axis xy equal; 
    colormap(ax1,'gray'); drawnow; hold on;
    set(ax1,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all')-500 max(mask_yrange,[],'all')+500]);
    %plot dummy lines to create a legend
    first_term = find(term_trace == 1,1,'first');
    % for l = 1:4
    %     pt(l) = plot(melmask.dated(first_term).x,melmask.dated(first_term).y,...
    %         '-','color',seas_cmap(l,:),'linewidth',1.5); hold on;
    % end

    %melange extent plotting: option 1 = plot the number of iceberg observations at each point to outline the
    %extent of the melange (instead of using melmask polygon)
    imx = single(im.x(min(xinds):max(xinds))); imy = single(im.y(min(yinds):max(yinds)));
    [Xgrid,Ygrid] = meshgrid(imx,imy); %convert the vector coordinates into matrices

    %loop through the DEMs & create a DEM mosaic
    disp('looping through the DEMs & loading elevations to a 3D matrix');
    melange_mats = dir([root_dir,site_abbrev,'/DEMs/*DEMfilled.mat']);
    for p = 1:length(melange_mats)
        DEM_name = melange_mats(p).name;
        load([root_dir,site_abbrev,'/DEMs/',DEM_name]);
        % figure; imagesc(M.DEM.x,M.DEM.y,M.DEM.z_filled); axis xy equal; hold on;
        % colormap jet;
        % plot(melmask.dated(p).x,melmask.dated(p).y,'-k','linewidth',2); hold on;
        [ZXgrid,ZYgrid] = meshgrid(M.DEM.x,M.DEM.y);

        %isolate the data in the melange mask
        DEMmask = zeros(size(ZXgrid));
        maskedz = M.DEM.z_filled.*M.mask.DEM;
        DEMmask(maskedz > zthresh) = 1;
        % figure; imagesc(M.DEM.x,M.DEM.y,DEMmask); axis xy equal; hold on;
        % colormap gray;

        %interpolate to the cropped Landsat image
        mel_z(:,:,p) = single(interp2(ZXgrid,ZYgrid,double(maskedz),Xgrid,Ygrid));
        mel_data(:,:,p) = single(interp2(ZXgrid,ZYgrid,double(DEMmask),Xgrid,Ygrid)); %interpolate elevations to the Landsat 8 coordinates
        clear Z* M DEMmask maskedz in;
    end
    melange_obs = sum(mel_data,3,"omitnan"); melange_z = (1026/(1026-917))*mean(mel_z,3,"omitnan"); 
    meltrans = melange_obs; meltrans(meltrans>0) = 1;
    z_cmap = cmocean('thermal',201); z_cmap(1,:) = [1 1 1];
    % % elev_cmap = colormap(pink(p)); 
    white_end = [231,212,232]/255; purple_end = [64,0,75]/255; colorgrad = (white_end-purple_end)/(p+1);
    for k = 1:p+1
        elev_cmap(k,:) = white_end-(k-1)*colorgrad;
    end
    % elev_cmap(1,:) = [1 1 1]; %white

    %plot the frequency of iceberg observations
    % figure1 = figure; set(gcf,'position',[50 50 1600 600]);
    ax2 = axes; 
    implot = imagesc(ax2,imx,imy,melange_z,'AlphaData', meltrans); axis xy equal; hold on;
    colormap(ax2,z_cmap); set(ax2,'color','none','visible','off'); set(gca,'clim',[0 200]);
    cbar = colorbar; cbar.YLabel.String = 'Average thickness (m)';
    % implot = imagesc(ax2,imx,imy,melange_obs/p,'AlphaData', meltrans); axis xy equal; hold on;
    % colormap(ax2,elev_cmap); set(ax2,'color','none','visible','off');
    % cbar = colorbar; cbar.YLabel.String = 'Iceberg cover frequency';
    cbar.FontSize = map_font-4; cbar.Location = "eastoutside";
    colormap(ax1,'gray');
    linkaxes([ax1 ax2]);
    plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',2);
    set(gca,'box','on','linewidth',1.5);

    %melange extent plotting
    for k = min(seaward_ext):max(inland_ext)+2
        %plot faint colored boxes for the DEM binning
        % fill([MP(j).Z.transectX(k,1),MP(j).Z.transectX(k,2),MP(j).Z.transectX(k+1,2),MP(j).Z.transectX(k+1,1),MP(j).Z.transectX(k,1)],...
        %     [MP(j).Z.transectY(k,1),MP(j).Z.transectY(k,2),MP(j).Z.transectY(k+1,2),MP(j).Z.transectY(k+1,1),MP(j).Z.transectY(k,1)],...
        %     tran_cmap(k-min(seaward_ext)+1,:),'FaceAlpha',0.5,'EdgeColor','none'); hold on;
        %plot box edges
        plot([MP(j).Z.transectX(k,1),MP(j).Z.transectX(k,2)],[MP(j).Z.transectY(k,1),MP(j).Z.transectY(k,2)],'--k','linewidth',2);
    end

    % %plot the terminus positions
    % termdists = MP(j).Z.termdist; termdists(term_trace == 0) = NaN;
    % [~,maxind] = max(termdists); [~,minind] = min(termdists);
    % plot(ax2,melmask.dated(minind).x,melmask.dated(minind).y,'-','color',[90,174,97]/255,'linewidth',3); hold on;
    % plot(ax2,melmask.dated(maxind).x,melmask.dated(maxind).y,'-','color',[90,174,97]/255,'linewidth',3); hold on;
    % disp('plotted min & max centerline terminus positions in green');
    % % for p = 1:size(melmask.dated,2)
    % %     if term_trace(p) == 1
    % %         plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
    % %             'k','linewidth',1.5); hold on;
    % %         % plot(melmask.dated(p).x,melmask.dated(p).y,'-','color',...
    % %         %     seas_cmap(zseas(p),:),'linewidth',1.5); hold on;
    % %         drawnow;
    % %     end
    % % end
    % % for k = 1:length(T)
    % %     if Tyrs(k) >= min(zyrs) && Tyrs(k) <= max(zyrs)
    % %     in = inpolygon(T(k).X,T(k).Y,melmask.uncropped.x,melmask.uncropped.y);
    % %     % plot(MP(j).T.X(in),MP(j).T.Y(in),'-','color',...
    % %     %     yr_cmap(Tyrs(k)-min(years)+1,:),'linewidth',1.5); hold on;
    % %     plot(T(k).X(in),T(k).Y(in),'-','color',...
    % %         seas_cmap(Tseas(k),:),'linewidth',1.5); hold on;
    % %     clear in; drawnow;
    % %     end
    % % end

    %add centerline & transect locations
    plot(ax2,MP(j).V.X(min(seaward_ext):end),MP(j).V.Y(min(seaward_ext):end),'-k','linewidth',2); hold on;
    if contains(big3,MP(j).name)
        plot(ax2,MP(j).V.X(min(seaward_ext):end),MP(j).V.Y(min(seaward_ext):end),'dk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
    else
        plot(ax2,MP(j).V.X(min(seaward_ext):end),MP(j).V.Y(min(seaward_ext):end),'sk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
    end
    set(ax2,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all')-500 max(mask_yrange,[],'all')+500]);
    % xlims = [min(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectX(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    % ylims = [min(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all'), max(MP(j).Z.transectY(min(seaward_ext):max(inland_ext)+1,:),[],'all')];
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    % set(ax1,'xlim',xlims,'ylim',ylims);
    xticks = get(ax1,'xtick'); yticks = get(ax1,'ytick');
    set(ax1,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',map_font);
    xlabel(ax1,'Easting (km)','fontsize',map_font); ylabel(ax1,'Northing (km)','fontsize',map_font);
    pos = get(gca,'position');
    if range(xlims) > 1.05*range(ylims) %short and fat map so plot the legend below
        % map_leg = legend(pt,season_names,'Location','southoutside',...
        %     'Orientation','horizontal');
        % leg_pos = get(map_leg,'position');
        % set(map_leg,'position',[leg_pos(1)+(0.5-mean([leg_pos(1) leg_pos(1)+leg_pos(3)])) 0.075 leg_pos(3) leg_pos(4)]);
        cbar.Location = "northoutside"; cbar.Orientation = 'horizontal';

        % if range(xlims)/range(ylims) < 1.5
            set(ax1,'position',[pos(1)+x_shift 0.2 0.9*pos(3) 0.9*pos(4)]); drawnow;
            set(ax2,'position',[pos(1)+x_shift 0.2 0.9*pos(3) 0.9*pos(4)]); drawnow;
        % else
        %     set(ax1,'position',[pos(1)+x_shift 0.2 pos(3) pos(4)]); drawnow;
        %     set(ax2,'position',[pos(1)+x_shift 0.2 pos(3) pos(4)]); drawnow;
        % end
    else %tall and thin map so plot the legend on the side
        % map_leg = legend(pt,season_names,'Location','eastoutside',...
        %     'Orientation','vertical');
        % leg_pos = get(map_leg,'position');
        % if leg_pos(1)+0.05 > pos(1)-0.05+pos(3)
        %     set(map_leg,'position',[leg_pos(1)+0.05 leg_pos(2) leg_pos(3) leg_pos(4)]);
        % else
        %     set(map_leg,'position',[pos(1)-0.05+pos(3) leg_pos(2) leg_pos(3) leg_pos(4)]);
        % end
        set(ax1,'position',[pos(1) pos(2) pos(3) pos(4)]); drawnow;
        set(ax2,'position',[pos(1) pos(2) pos(3) pos(4)]); drawnow;
    end
    saveas(map_fig,[root_dir,site_abbrev,'/',site_abbrev,'-',suffix,'.png'],'png'); %save the image
    exportgraphics(map_fig,[root_dir,site_abbrev,'/',site_abbrev,'-',suffix,'.tif'],Resolution=600);
    close(map_fig);
    % uiwait %advance only after figure is closed

    %for Alison Glacier for 20110610, create a two subpanel figure to
    %demonstrate how packing density and the bergy bit misfit are
    %calculated with the map on the top & the size distribution on the bottom
    if contains(site_abbrev,'ASG')
        ex_fig = figure; set(ex_fig,'position',[50 50 900 1200]);

        %reproduce the site overview map
        subo = subplot(3,1,1); 
        %plot the image
        imagesc(subo,im.x(min(xinds):max(xinds)),im.y(min(yinds):max(yinds)),im_subset); axis xy equal;
        colormap(subo,'gray'); drawnow; hold on;
        set(subo,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all') max(mask_yrange,[],'all')-1000]);
        subo_pos = get(gca,'position'); set(subo,'position',[subo_pos(1) subo_pos(2)-0.1 subo_pos(3) subo_pos(4)]);
        %add the average thickness
        ax2 = axes;
        implot = imagesc(ax2,imx,imy,melange_z,'AlphaData', meltrans); axis xy equal; hold on;
        colormap(ax2,z_cmap); set(ax2,'color','none','visible','off'); set(gca,'clim',[0 200]);
        cbar = colorbar; cbar.YLabel.String = 'Average thickness (m)';
        cbar.FontSize = 16; cbar.Location = "eastoutside";
        colormap(subo,'gray');
        linkaxes([subo ax2]); set(gca,'box','on','linewidth',1.5);
        plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',2);
        %melange extent plotting
        for k = min(seaward_ext):max(inland_ext)+2
            plot([MP(j).Z.transectX(k,1),MP(j).Z.transectX(k,2)],[MP(j).Z.transectY(k,1),MP(j).Z.transectY(k,2)],'--k','linewidth',2);
        end
        %add centerline & transect locations
        plot(ax2,MP(j).V.X(min(seaward_ext):end),MP(j).V.Y(min(seaward_ext):end),'-k','linewidth',2); hold on;
        if contains(big3,MP(j).name)
            plot(ax2,MP(j).V.X(min(seaward_ext):end),MP(j).V.Y(min(seaward_ext):end),'dk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
        else
            plot(ax2,MP(j).V.X(min(seaward_ext):end),MP(j).V.Y(min(seaward_ext):end),'sk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
        end
        set(ax2,'xlim',[min(mask_xrange,[],'all')-500 max(mask_xrange,[],'all')+500],'ylim',[min(mask_yrange,[],'all') max(mask_yrange,[],'all')-1000]);
        xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
        xticks = get(subo,'xtick'); yticks = get(subo,'ytick');
        set(subo,'xticklabels',[],'fontsize',16); %xlabel(subo,'Easting (km)','fontsize',16); 
        set(subo,'yticklabels',yticks/1000,'fontsize',16); ylabel(subo,'Northing (km)','fontsize',16); 
        if range(xlims) > 1.05*range(ylims) %short and fat map so plot the legend below
            cbar.Location = "eastoutside"; cbar.Orientation = 'vertical';
            set(subo,'position',[0 0.7 1 (subo_pos(4)/subo_pos(3))]); drawnow;
            set(ax2,'position',[0 0.7 1 (subo_pos(4)/subo_pos(3))]); drawnow;
        else %tall and thin map so plot the legend on the side
            set(subo,'position',[subo_pos(1) subo_pos(2) subo_pos(3) subo_pos(4)]); drawnow;
            set(ax2,'position',[subo_pos(1) subo_pos(2) subo_pos(3) subo_pos(4)]); drawnow;
        end
        subo_pos = get(subo,'position');

        %load the DEM of interest for the example map
        melange_mats = dir([root_dir,site_abbrev,'/DEMs/*DEMfilled.mat']);
        for p = 1:length(melange_mats)
            DEM_name = melange_mats(p).name;
            if contains(DEM_name,'20110610')
                load([root_dir,site_abbrev,'/DEMs/',DEM_name]);
                [ZXgrid,ZYgrid] = meshgrid(M.DEM.x,M.DEM.y);

                %create the packing density mask
                packmask = zeros(size(ZXgrid));
                maskedz = M.DEM.z_filled.*M.mask.DEM;
                packmask(maskedz <= zthresh) = 1;

                clear Z*;
            end
        end

        %plot the map
        subm = subplot(3,1,2);
        maskedz(isnan(maskedz)) = -0.1; axa = gca;
        imagesc(subm,M.DEM.x,M.DEM.y,(1026/(1026-917))*maskedz); axis xy equal;
        % berg_cmap = colormap(gray(501)); berg_cmap(1,:) = [1 1 1];
        % colormap(subm,berg_cmap); drawnow; hold on;
        colormap(subm,z_cmap); set(subm,'clim',[0 200]); drawnow; hold on;
        % set(subm,'xlim',[-332500 -320500],'ylim',[-1646500 -1638500]);
        % berg_cbar = colorbar; %berg_cbar.Ticks = [0:50:200]; 
        % berg_cbar.Label.String = 'Thickness (m)';
        % berg_cbar.FontSize = 16; berg_cbar.Location = 'northoutside';
        %format the map
        % xticks = get(subm,'xtick'); yticks = get(subm,'ytick');
        % set(subm,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',14);
        % xlabel(subm,'Easting (km)','fontsize',14); ylabel(subm,'Northing (km)','fontsize',14);
        set(subm,'xlim',xlims,'ylim',ylims); set(gca,'box','on','linewidth',1.5);
        xticks = get(subm,'xtick'); yticks = get(subm,'ytick');
        set(subm,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
        xlabel(subm,'Easting (km)','fontsize',16); ylabel(subm,'Northing (km)','fontsize',16);
        subm_pos = get(subm,'position');
        set(subm,'position',[0 0.4 1 (subo_pos(4)/subo_pos(3))]);

        %load the size distribution for near the terminus for the date
        D = readtable([root_dir,site_abbrev,'/',site_abbrev,'-20110610-iceberg-distribution-subsets.csv'],"VariableNamingRule","preserve");
        v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes = sum(~isnan(berg_nos),1);
        inland_idx = find(size_classes>0,1,'last');
        n1 = berg_nos(:,inland_idx)./dv1;
        clear D;
        %filter the data & fit the curve
        v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
        v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs
        %enter unfiltered data into model (filtering happens there)
        tempfig = figure; drawnow;
        figure(tempfig);
        [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
        n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model
        close(tempfig);

        %add the mask for bergy bits
        packmask(maskedz==-0.1) = -1; packmask(maskedz==0) = -1;
        packalpha = 0.8*packmask;
        figure(ex_fig); subplot(subm);
        axb = axes; %mask_cmap = [1,1,1; 0,0,0; 251/255,180/255,185/255];
        mask_cmap = [1,1,1; 0,0,0; 0.85 0.85 0.85];
        implot = imagesc(axb,M.DEM.x,M.DEM.y,packmask,'AlphaData', packalpha); axis xy equal; hold on;
        colormap(axb,mask_cmap); set(axb,'color','none','visible','off');
        % set(axb,'xlim',[-332500 -320500],'ylim',[-1646500 -1638500]);
        set(axb,'xlim',xlims,'ylim',ylims);
        % colormap(subm,berg_cmap);
        colormap(subm,z_cmap);
        map_pos = get(subm,'position'); set(subm,'position',[0 0.4 1 (subo_pos(4)/subo_pos(3))]);
        map_pos = get(subm,'position'); set(axb,'position',map_pos);
        clear packmask maskedz M;
        % %add the centerline
        % plot(axb,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'-k','markersize',6,'linewidth',2); hold on;
        % plot(axb,MP(j).V.X(min(seaward_ext):max(inland_ext)+1),MP(j).V.Y(min(seaward_ext):max(inland_ext)+1),'sk','markersize',6,'linewidth',1,'markerfacecolor','k'); hold on;
        %load the transects
        S = shaperead([root_dir,site_abbrev,'/shapefiles/',site_abbrev,'-transects_1000m.shp']);
        clear *_ext;

        %plot the transect bounds for the subset
        % plot(melmask.uncropped.x,melmask.uncropped.y,'-k'); hold on;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx).X(1)).^2 + (melmask.uncropped.y-S(inland_idx).Y(1)).^2);
        [~,idx] = min(end_dist); end_idx(1) = idx; clear end_dist idx;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx).X(end-1)).^2 + (melmask.uncropped.y-S(inland_idx).Y(end-1)).^2);
        [~,idx] = min(end_dist); end_idx(2) = idx; clear end_dist idx;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx+1).X(1)).^2 + (melmask.uncropped.y-S(inland_idx+1).Y(1)).^2);
        [~,idx] = min(end_dist); start_idx(1) = idx; clear end_dist idx;
        end_dist = sqrt((melmask.uncropped.x-S(inland_idx+1).X(end-1)).^2 + (melmask.uncropped.y-S(inland_idx+1).Y(end-1)).^2);
        [~,idx] = min(end_dist); start_idx(2) = idx; clear end_dist idx;
        plot([S(inland_idx).X(1:end-1),melmask.uncropped.x(end_idx(2)+1:1:start_idx(2))',S(inland_idx+1).X(end-1:-1:1),melmask.uncropped.x(start_idx(1)+1:1:end_idx(1)-1)',S(inland_idx).X(1)],...
            [S(inland_idx).Y(1:end-1),melmask.uncropped.y(end_idx(2)+1:1:start_idx(2))',S(inland_idx+1).Y(end-1:-1:1),melmask.uncropped.y(start_idx(1)+1:1:end_idx(1)-1)',S(inland_idx).Y(1)],'-','linewidth',2,...
            'color',seas_cmap(2,:)); hold on;
        plot(melmask.uncropped.x,melmask.uncropped.y,'-k','linewidth',2);
        clear S *_idx;

        %shift the colorbar so it aligns with both maps
        cbar_pos = cbar.Position; cbar.Position = [cbar_pos(1)+0.02 0.5 cbar_pos(3) cbar_pos(4)+0.1];

        %plot the data
        figure(ex_fig);
        subp = subplot(3,1,3);
        plot(subp,log10([v1(zthresh+1),v1(zthresh+1)]),log10([10^-12 10^3]),'-.','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
        plot(subp,log10(v1),log10(n1),'-','LineWidth',2,'Color',seas_cmap(2,:)); hold on;
        plot(subp,log10(v1),log10(n1),'s','LineWidth',1,'Color',seas_cmap(2,:),'markersize',5,'markerfacecolor',seas_cmap(2,:)); hold on;
        plot(subp,log10(v1),log10(n1_mod),'--','LineWidth',2,'Color',seas_cmap(2,:)); hold on;
        %format the plot
        subp.FontSize = 16;
        subp.YLim = [-5 2]; subp.YTick = [-4:2:2]; %axl.YLim = [10^-10 10^3];
        subp.XLim = [1.5 5]; subp.XTick = [2, 3, 4, 5];
        subp.XTickLabel = compose('%.0e',[10^2, 10^3, 10^4, 10^5]);
        subp.YTickLabel = compose('%.0e', [10^-4, 10^-2, 0, 10^2]);
        set(gca,'box','on','linewidth',1.5);
        xlabel('Surface area (m^2)','fontsize',16);
        ylabel('Normalized iceberg count','fontsize',16);
        plot_pos = get(gca,'position');
        set(subp,'position',[0.165 0.10 0.67 0.2]);

        %save the figure
        saveas(ex_fig,[root_dir,site_abbrev,'/',site_abbrev,'-bergy-bit_example-subplots.png'],'png'); %save the image
        exportgraphics(ex_fig,[root_dir,site_abbrev,'/',site_abbrev,'-bergy-bit_example-subplots.tif'],Resolution=600);
        close(ex_fig);
        clear inland_idx v n dv n1 v1 dv1;
    end

    clear imx imy Xgrid Ygrid DEMmask;
    clear xlims ylims pt elev_cmap maskedz meltrans mel_data mel_z melange_* termdists maxind minind mask_*range;
    clear berg_* im* LCdir *_ext leg* map_fig map_leg melmask zmos on pos ref* seaward_ext size_classes sort_ind T T_inds term_* Tyrs Tmos *lims zyrs zdate zseas Tseas tran_cmap colorgrad term;
end
close all;
disp('Done creating site maps');
end


%% plot iceberg size distributions as a gif & as subplots
%distribution for a site, then adds the fall distribution (if it exists)
%with arrows showing change between seasons for small (~100 m^2) and big
%(~10^5 m^2) icebergs

%find size indices for arrows
si = find(MP(1).D.area <= 100,1,'last');
bi = find(MP(1).D.area <= 10^5,1,'last');

if zcutoff == 3
%create the figure
dist_fig = figure; set(dist_fig,'position',[50 850 800 600]);
loglog(MP(1).D.area,MP(geo_ind(1)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',16); grid on; 
xlabel('Surface area (m^2)','fontsize',16); ylabel('Normalized iceberg count','fontsize',16); 
drawnow;
nimages = 1;
for j = 1:length(geo_ind)
    if sum(~isnan(MP(geo_ind(j)).D.bergs(2,:))) > 0
        %plot spring
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
        title(MP(geo_ind(j)).name);
        drawnow;
        frame = getframe(dist_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
        
        for p = 3:4
            %plot summer and/or spring then fall depending on data
            if sum(~isnan(MP(geo_ind(j)).D.bergs(p,:))) > 0
                loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
                title(MP(geo_ind(j)).name);
                drawnow;
                
                %later distribution
                loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(p,:),'-','color',seas_cmap(p,:),'linewidth',2); hold on;
                drawnow;
                %             frame = getframe(dist_fig);
                %             gif_im{nimages} = frame2im(frame); nimages = nimages+1;
                
                %add arrows for small iceberg abundance change
                sp1 = [MP(1).D.area(si) MP(geo_ind(j)).D.bergs(2,si)];
                sp2 = [MP(1).D.area(si) MP(geo_ind(j)).D.bergs(p,si)];
                sdp = sp2-sp1;
                quiver(sp1(1),sp1(2),sdp(1),sdp(2),0,'linewidth',2,'color','k'); hold on;
                if sdp(2) > 0
                    plot(sp2(1),sp2(2)-0.2*sp2(2),'^k','markerfacecolor','k'); hold on;
                else
                    plot(sp2(1),sp2(2)+0.2*sp2(2),'vk','markerfacecolor','k'); hold on;
                end
                clear sp1 sp2 dsp;
                %add arrows for large iceberg abundance change
                sp1 = [MP(1).D.area(bi) MP(geo_ind(j)).D.bergs(2,bi)];
                sp2 = [MP(1).D.area(bi) MP(geo_ind(j)).D.bergs(p,bi)];
                sdp = sp2-sp1;
                quiver(sp1(1),sp1(2),sdp(1),sdp(2),0,'linewidth',2,'color','k'); hold on;
                if sdp(2) > 0
                    plot(sp2(1),sp2(2)-0.2*sp2(2),'^k','markerfacecolor','k'); hold on;
                else
                    plot(sp2(1),sp2(2)+0.2*sp2(2),'vk','markerfacecolor','k'); hold on;
                end
                clear sp1 sp2 dsp;
                drawnow;
                frame = getframe(dist_fig);
                gif_im{nimages} = frame2im(frame); nimages = nimages+1;
                
                %clear the plot
                cla;
            end
        end
        cla;
        title('');
        frame = getframe(dist_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    end
end
close;

filename = [root_dir,'GrIS-melange_seasonal-iceberg-size-distributions_timeseries.gif']; % Specify the output file name
for idx = 1:nimages-1
    [A,map] = rgb2ind(gif_im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif",LoopCount=Inf, ...
                DelayTime=1)
    else
        imwrite(A,map,filename,"gif",WriteMode="append", ...
                DelayTime=1)
    end
end
close all; clear A; clear frame gif_im map;

%plot all the seasonal distributions on a single figure
dist_fig = figure; set(dist_fig,'position',[50 50 1800 500]);
sub1 = subplot(1,3,1); sub2 = subplot(1,3,2); sub3 = subplot(1,3,3);
for j = 1:length(geo_ind)
        %plot spring
        subplot(sub1);
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(2,:),'-','color',seas_cmap(2,:),'linewidth',2); hold on;
        drawnow;
        
        %plot summer
        subplot(sub2);
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(3,:),'-','color',seas_cmap(3,:),'linewidth',2); hold on;
        
        %plot fall
        subplot(sub3);
        loglog(MP(geo_ind(j)).D.area,MP(geo_ind(j)).D.bergs(4,:),'-','color',seas_cmap(4,:),'linewidth',2); hold on;
        drawnow;

end
subplot(sub1);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'fontsize',16); grid on; 
text(10^5.8,10^-1,'spring','fontsize',16)
xlabel('Surface area (m^2)','fontsize',16); ylabel('Normalized iceberg count','fontsize',16);
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub2);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',16); grid on;
text(10^5.8,10^-1,'summer','fontsize',16)
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
subplot(sub3);
set(gca,'ylim',[10^-12 1],'xlim',[10^1 10^7],'yticklabel',[],'fontsize',16); grid on; 
text(10^5.8,10^-1,'fall','fontsize',16)
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.03 1.2*pos(3) pos(4)]);
drawnow;
saveas(dist_fig,[root_dir,'GrIS-melange_iceberg-distribution_loglog.png'],'png'); %save the plot
close all;
end

%% plot seasonal size distributions with empirical fits to demonstrate bergy bit misfits

%specify if you want to zoom in the plots to focus on the power-law parts
%of the size distributions
suffix = 'plzoom_loglog'; %'plzoom_loglog' to focus on iceberg areas <10^5m^2, 'fulldist_loglog' for full

% %Size distribution parameters (only needed if you didn't run the sectionabove with the same info)
% nthresh = 1e-6; % set small number bin cutoff (n1 must be greater than this value)
% zthresh = 3; %set small size bin cutoff (freeboard must exceed this value)
% rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)
% ARcomp.best.autoALL = 2; % iceberg aspect ratio 
% vthresh = (1/4)*pi*((rho_sw/(rho_sw-rho_i))*ARcomp.best.autoALL.*zthresh).^2; %don't include this bin size or smaller in curve fitting
% dplawthresh = 10^5; % upper bound on the intercept for the dummy powerlaw
% norm_type = 2; % toggle between L2, max, and log norm using 2, Inf, and 'log'
% normalize_exp = 1.5; % Increase to weight residuals towards end of the curve, minimum = 1

%create a seasonal profile figure with a top row containing profiles near the terminus and a bottom row containing seaward profiles
subdist_fig = figure; set(subdist_fig,'position',[50 50 1200 800]);
%subplots for inland sample location
subia = subplot(2,3,1); subib = subplot(2,3,2); subic = subplot(2,3,3); 
%subplots for seaward sample location
subsa = subplot(2,3,4); subsb = subplot(2,3,5); subsc = subplot(2,3,6); 
%letters for the plots
alphabet = ['a','b','c','d','e','f'];
for j = 1:length(geo_ind)
    seasons = MP(geo_ind(j)).D.months; site_abbrev = MP(geo_ind(j)).name;
    % disp(sitenames(geo_ind(j),:)); 
    
    %isolate the dates for data sorting
    for p = 1:length(MP(geo_ind(j)).Z.date)
        term_trace(p) = MP(geo_ind(j)).Z.termflag(p);
        % zdate(p) = convert_to_decimaldate(char(MP(geo_ind(j)).Z.date(p)));
        % datest = datetime(MP(geo_ind(j)).Z.date{p},'InputFormat','yyyyMMdd');
        % zyrs(p) = year(datest); zmos(p) = month(datest);
        % clear datest;
    end

    %loop through subset size distributions & find the first and last
    %non-NaN columns to identify melange extent
    Dsubs = dir([root_dir,MP(geo_ind(j)).name,'/',MP(geo_ind(j)).name,'*-iceberg-distribution-subsets.csv']);
    size_classes = [];
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(geo_ind(j)).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");
        name_split = split(Dsubs(p).name,'-',2);
        berg_dates(p) = datetime(name_split{2},'Inputformat','yyyyMMdd'); clear name_split;
        bmos(p) = month(berg_dates(p)); byrs(p) = year(berg_dates(p));

        %identify observational limits along the centerline
        berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
        size_classes(p,:) = sum(~isnan(berg_nos),1);
        % seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
        inland_ext(p) = find(size_classes(p,:)>0,1,'last')+1; inland_idx(p) = inland_ext(p);
        clear berg_nos;

        %fill-in inland boundary reference based on terminus observations
        %from images if DEM did not overlap the terminus
        if term_trace(p) == 0
            if ~isnan(MP(geo_ind(j)).Z.termdist(p))
                dists = sqrt((MP(geo_ind(j)).Z.termX(p)-MP(geo_ind(j)).V.X).^2 + (MP(geo_ind(j)).Z.termY(p)-MP(geo_ind(j)).V.Y).^2);
                [~,ia] = sort(dists);
                inland_idx(p) = max(ia(1:2)); term_trace(p) = 1; %use terminus traces from imagery
                clear dists ia;
            else
                inland_idx(p) = NaN;
                disp(['Need terminus data for ',MP(geo_ind(j)).name,' ',char(MP(geo_ind(j)).Z.date(p))])
            end
        end
        
    end
    clear size_classes;

    %now analyze the size distributions
    for p = 1:length(Dsubs)
        D = readtable([root_dir,MP(geo_ind(j)).name,'/',Dsubs(p).name],"VariableNamingRule","preserve");

        %extract size distribution profiles
        if p == 1 %initialize the data compilation
            v1 = table2array(D(:,1)); dv1 = table2array(D(:,2));
        end
        %first 2 columns are area & bin width, so 3+ are data from points along the centerline
        n1 = NaN(length(v1),max(inland_idx)-1);
        if ~isnan(inland_idx(p))
            berg_nos = fliplr(table2array(D(:,3:(inland_idx(p)-1)+2)));
            berg_nos(berg_nos==0) = NaN;
            n1(:,1:size(berg_nos,2)) = berg_nos./dv1;

            size_classes = sum(~isnan(berg_nos),1);
            seaward_ext(p) = find(size_classes>0,1,'last'); %seaward limit index in moving reference frame relative to the terminus
        end
        clear D;
        %stack distributions in a structure
        binned_n1_stack(:,:,p) = n1;
        clear n1 berg_nos size_classes; 
    end
    %identify seaward sample location
    seaward_meanidx = round(mean(seaward_ext));


    %pasted from above...
    %create seasonal average distributions for each fjord bin
    n1wi_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    n1sp_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    n1su_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    n1fa_stack = NaN(size(binned_n1_stack,1),size(binned_n1_stack,2),length(years));
    for p = 1:length(years)
        yr_idx = find(byrs == years(p));
        if ~isempty(yr_idx)

            %isolate the distributions for that year
            for k = 1:length(yr_idx)
                n1_temp(:,:,k) = binned_n1_stack(:,:,yr_idx(k));
            end

            %calculate seasonal averages across the year
            mos_yr = bmos(yr_idx);
            if ~isempty(ismember(mos_yr,seasons(1,:))) %winter
                n1wi_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(1,:))==1),3);
            end
            if ~isempty(ismember(mos_yr,seasons(2,:))) %spring
                n1sp_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(2,:))==1),3);
            end
            if ~isempty(ismember(mos_yr,seasons(3,:))) %summer
                n1su_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(3,:))==1),3);
            end
            if ~isempty(ismember(mos_yr,seasons(4,:))) %fall
                n1fa_stack(:,:,p) = nanmean(n1_temp(:,:,ismember(mos_yr,seasons(4,:))==1),3);
            end

            clear mos_yr n1_temp;
        end
        clear yr_idx
    end
    %calculate seasonal averages across all years & fit curves
    n1_seasmean(:,:,1) = nanmean(n1wi_stack,3);
    n1_seasmax(:,:,1) = (nanmean(n1wi_stack,3)+std(n1wi_stack,0,3,'omitnan')); n1_seasmin(:,:,1) = (nanmean(n1wi_stack,3)-std(n1wi_stack,0,3,'omitnan'));
    n1_seasmean(:,:,2) = nanmean(n1sp_stack,3);
    n1_seasmax(:,:,2) = (nanmean(n1sp_stack,3)+std(n1sp_stack,0,3,'omitnan')); n1_seasmin(:,:,2) = (nanmean(n1sp_stack,3)-std(n1sp_stack,0,3,'omitnan'));
    n1_seasmean(:,:,3) = nanmean(n1su_stack,3);
    n1_seasmax(:,:,3) = (nanmean(n1su_stack,3)+std(n1su_stack,0,3,'omitnan')); n1_seasmin(:,:,3) = (nanmean(n1su_stack,3)-std(n1su_stack,0,3,'omitnan'));
    n1_seasmean(:,:,4) = nanmean(n1fa_stack,3);
    n1_seasmax(:,:,4) = (nanmean(n1fa_stack,3)+std(n1fa_stack,0,3,'omitnan')); n1_seasmin(:,:,4) = (nanmean(n1fa_stack,3)-std(n1fa_stack,0,3,'omitnan'));
    clear n1wi_stack n1sp_stack n1su_stack n1fa_stack;

    %create seasonal average distributions for each fjord bin
    for k = 1:4
        for l = [1,seaward_meanidx]
            n1 = n1_seasmean(:,l,k);
            if l == 1
                n1_inlandseas(:,j,k) = n1;
            else
                n1_seawardseas(:,j,k) = n1;
            end
            clear n1;
        end
    end
    %concatenate the profiles for each season so that they can be plotted
    %as filled polygons instead of individual lines
    %seaward
    n1wi(:,j,2) = n1_seawardseas(:,j,1); n1sp(:,j,2) = n1_seawardseas(:,j,2); 
    n1su(:,j,2) = n1_seawardseas(:,j,3); n1fa(:,j,2) = n1_seawardseas(:,j,4);
    %inland
    n1wi(:,j,1) = n1_inlandseas(:,j,1); n1sp(:,j,1) = n1_inlandseas(:,j,2); 
    n1su(:,j,1) = n1_inlandseas(:,j,3); n1fa(:,j,1) = n1_inlandseas(:,j,4);

    clear berg_dates bmos byrs berg_nos size_classes seaward_* inland_* binned_n1_stack n1_seasmean term_trace;
    clear n1_seasm*;
end
log_v = log10(v1); 

%apply the EBC fitc function to the average of all seasonal profiles in each
%location & plot the result as a dashed line
sample_locs = {'terminus','seaward'};
for l = 1:2
    if l == 1
        disp('plotting profiles for near-terminus sample location')
    else
        disp('plotting profiles for seaward sample location')
    end

    %2 = spring, 3 = summer, 4 = fall
    for k = 2:4
        if k <= 2
            fill_alpha = 0.2;
        else
            fill_alpha = 0.4;
        end

        if l == 1
            n1 = nanmedian(n1_inlandseas(:,:,k),2);
        else
            n1 = nanmedian(n1_seawardseas(:,:,k),2);
        end

        %filter to decide if there are enough points to fit lines
        v = v1(n1>nthresh); n = n1(n1>nthresh); dv = dv1(n1>nthresh); % TOGGLE TO REMOVE SMALLEST COUNTS
        n = n(v>vthresh); dv = dv(v>vthresh); v = v(v>vthresh); %remove smallest size classes
        v = v(~isnan(n)); dv = dv(~isnan(n)); n = n(~isnan(n)); % Remove NaNs

        % Fit E-BC fragmentation model to the data and plot
        if sum(~isnan(v)) > 12
            %enter unfiltered data into model (filtering happens there)
            tempfig = figure;
            [alpha,c1,c2,c3,c4,~,~] = EBC_fragmentation_curve(v1, n1, dv1, norm_type, vthresh, nthresh, 1, normalize_exp); % fit Eq. 1
            n1_mod = EBC_model([c1,c2,alpha,c3,c4],v1); % grab model
            close(tempfig);

            %plot the size distribution fit
            figure(subdist_fig);
            if l == 1
                if k == 2
                    subplot(subia); log_n = log10(n1sp); 
                elseif k == 3
                    subplot(subib); log_n = log10(n1su);
                else
                    subplot(subic); log_n = log10(n1fa);
                end
                yyaxis left; axl = gca; pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)-0.03 1.075*pos(3) 1.05*pos(4)]);

                %plot dummy fill polygons for the legend
                if l == 1 && k == 2
                    pf(1) = fill(gca,[10,50,50,10,10],[0.9e2,0.9e2,1.1e2,1.1e2,0.9e2],...
                        seas_cmap(2,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on;
                    pf(2) = fill(gca,[10,50,50,10,10],[0.9e2,0.9e2,1.1e2,1.1e2,0.9e2],...
                        seas_cmap(3,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on;
                    pf(3) = fill(gca,[10,50,50,10,10],[0.9e2,0.9e2,1.1e2,1.1e2,0.9e2],...
                        seas_cmap(4,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on;
                end

                %plot the data
                plot(log10([v1(zthresh+1),v1(zthresh+1)]),log10([10^-12 10^3]),'-.','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                plot(log10(v1),log10(n1),'.-','LineWidth',2,'Color','k'); hold on;
                plot(log10(v1),log10(n1_mod),'--','LineWidth',2,'Color','k'); hold on;
                % loglog([v1(zthresh+1),v1(zthresh+1)],[10^-12 10^3],'-','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                % loglog(v1,n1,'.-','LineWidth',2,'Color','k'); hold on;
                % loglog(v1,n1_mod,'--','LineWidth',2,'Color','k'); hold on;

                %format the subplot
                axl.YAxis(1).Color = 'k'; set(axl,'box','on'); axl.LineWidth = 1.5;
                if contains(suffix,'zoom')
                    axl.YLim = [-6 3]; axl.YTick = [-6:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 5]; axl.XTick = [2, 3, 4, 5];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^3, 10^4, 10^5]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(3.5,2.3,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',18);
                else
                    axl.YLim = [-8 3]; axl.YTick = [-8:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 6.5]; axl.XTick = [2, 4, 6];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^4, 10^6]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-8, 10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(4.2,2.2,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',18);
                end
                axl.FontSize(1) = 18; %grid on;
                %axl.XLim = [10^1 10^7]; axl.XTick = [10^2, 10^4, 10^6];
                % set(gca,'ylim',[10^-10 10^3],'xlim',[10^1 10^7],'xtick',[10^2, 10^4, 10^6],'fontsize',18); grid on;
                % text(10^5,10^2,[alphabet(k*l-1),') ',char(season_names(k))],'fontsize',18)
                if k == 2
                    ylabel('Normalized iceberg count','fontsize',18);
                end
                
            else
                if k == 2
                    subplot(subsa); log_n = log10(n1sp); 
                elseif k == 3
                    subplot(subsb); log_n = log10(n1su);
                else
                    subplot(subsc); log_n = log10(n1fa);
                end
                yyaxis left; axl = gca; pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.075*pos(3) 1.05*pos(4)]);

                %plot the data
                plot(log10([v1(zthresh+1),v1(zthresh+1)]),log10([10^-12 10^3]),'-.','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                plot(log10(v1),log10(n1),'.-','LineWidth',2,'Color','k'); hold on;
                plot(log10(v1),log10(n1_mod),'--','LineWidth',2,'Color','k'); hold on;
                % loglog([v1(zthresh+1),v1(zthresh+1)],[10^-12 10^3],'-','LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
                % loglog(v1,n1,'.-','LineWidth',2,'Color','k'); hold on;
                % loglog(v1,n1_mod,'--','LineWidth',2,'Color','k'); hold on;

                %format the subplot
                axl.YAxis(1).Color = 'k'; set(axl,'box','on'); axl.LineWidth = 1.5;
                if contains(suffix,'zoom')
                    axl.YLim = [-6 3]; axl.YTick = [-6:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 5]; axl.XTick = [2, 3, 4, 5];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^3, 10^4, 10^5]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(3.5,2.3,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',18);
                else
                    axl.YAxis(1).Color = 'k'; set(axl,'box','on'); axl.LineWidth = 1.5;
                    axl.YLim = [-8 3]; axl.YTick = [-8:2:2]; %axl.YLim = [10^-10 10^3];
                    axl.XLim = [1.5 6.5]; axl.XTick = [2, 4, 6]; %axl.XLim = [10^1 10^7]; axl.XTick = [10^2, 10^4, 10^6];
                    axl.XTickLabel = compose('%.0e',[10^2, 10^4, 10^6]);
                    if k == 2
                        axl.YTickLabel = compose('%.0e', [10^-8, 10^-6, 10^-4, 10^-2, 0, 10^2]);
                    else
                        axl.YTickLabel = [];
                    end
                    text(4.2,2.2,[alphabet((l-1)*3 + (k-1)),') ',char(sample_locs(l))],'fontsize',18);
                end
                axl.FontSize(1) = 18; %grid on;
                % set(gca,'ylim',[10^-10 10^3],'xlim',[10^1 10^7],'xtick',[10^2, 10^4, 10^6],'fontsize',18); 
                % text(10^5,10^2,[alphabet((k-1)+3),') ',char(season_names(k))],'fontsize',18)
                xlabel('Surface area (m^2)','fontsize',18);
                if k == 2
                    ylabel('Normalized iceberg count','fontsize',18);
                end
            end
            %add the raw data
            yyaxis right; axr = gca; axr.XScale = 'linear';
            pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
            maxn = max(log_n(:,:,l),[],2); minn = min(log_n(:,:,l),[],2);
            fill(gca,[log_v(~isnan(maxn)); flipud(log_v(~isnan(minn)))],[maxn(~isnan(maxn)); flipud(minn(~isnan(minn)))],...
                seas_cmap(k,:),'FaceAlpha',fill_alpha,'EdgeColor','none'); hold on;
            if contains(suffix,'zoom')
                set(axr,'ylim',[-6 3],'xlim',[1.5 5],'yticklabel',[]);
            else
                set(axr,'ylim',[-8 3],'xlim',[1.5 6.5],'yticklabel',[]);
            end
            clear maxn minn;
            axl.YAxis(2).Color = 'k';

            clear n1_mod v dv n n1 alpha c1 c2 c3 c4;
        end
        drawnow;

        %add the legend
        if l == 1 && k == 2
            seas_leg = legend(pf,season_names(2:4));
            set(seas_leg,'location','southoutside','orientation','horizontal')
            set(gca,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
            set(axl,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
            set(axr,'position',[pos(1) pos(2) 1.05*pos(3) pos(4)]);
        end
        clear axl axr;
    end
    disp('moving on')
end
set(seas_leg,'position',[0.43 0.93 0.2,0.03]);
saveas(subdist_fig,[root_dir,'GrIS-melange_iceberg-distribution-variability_',num2str(zcutoff),'m-zcutoff_',suffix,'.png'],'png'); %save the plot
exportgraphics(subdist_fig,[root_dir,'GrIS-melange_iceberg-distribution-variability_',num2str(zcutoff),'m-zcutoff_',suffix,'.tif'],Resolution=600);

