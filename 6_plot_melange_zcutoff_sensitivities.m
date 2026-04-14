%%% Create results of GrIS melange sensitivity tests
clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/ArcticMappingTools/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/miscellaneous/melange-characterization-code/');

%specify root path
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';

%customize visualization
years = 2011:1:2023; yr_cmap = cmocean('matter',length(years)+1); yr_cmap = yr_cmap(2:end,:);
seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MAM','JJA','SON'};
seas_cmap = [5,113,176; 202,0,32; 244,165,130; 146,197,222]/255; %blue/orange/red
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


%% display sensitivity test results figure & the seasonal information for each site (GRL paper Table)
disp('Seasonal statistics for all sites:'); close all; drawnow;

%locate all the results for different zcutoff values
MPfiles = dir('GrIS-melange-characteristics_*m-zcutoff.mat'); 

%create a colormap for the sensitivity test results
zcut_cmap = [194,165,207; 118,42,131; 0,0,0; 27,120,55; 166,219,160]./255;

%grab the data for each zctuoff, plot, & aggregate for stats
zcutfig = figure; set(zcutfig,'position',[50 850 1200 900]);
%sensitivity test results
subH = subplot(2,2,1); subB = subplot(2,2,2); 
ax1 = gca; ax2 = axes; set(ax2,'position',get(ax1,'position')); 
%best results color-coded by season
subHs = subplot(2,2,3); subBs = subplot(2,2,4); 
ax1s = gca; ax2s = axes; set(ax2s,'position',get(ax1s,'position')); 
for i = 1:length(MPfiles)
    load([root_dir,MPfiles(i).name]);
    disp(MPfiles(i).name);
    for j = 1:length(geo_ind)
        %grab date info
        for p = 1:length(MP(geo_ind(j)).Z.date)
            zdate(p) = convert_to_decimaldate(char(MP(geo_ind(j)).Z.date(p)));
            datest(p,:) = datetime(MP(geo_ind(j)).Z.date{p},'InputFormat','yyyyMMdd');
            dateout(p,:) = datestr(datest(p,:),'yyyy-mm-dd');
            zyrs(p) = year(datest(p,:)); zmos(p) = month(datest(p,:));
            % clear datest;
        end

        %site abbreviation and name
        % disp(MP(geo_ind(j)).name);
        disp(char(geo_names(j)));

        %seasonal DEM dates, near-terminus thickness, buttressing
        for k = 1:4
            if i == 3
                disp(['  ',char(season_names(k))]);

                %DEM dates
                seas_refs = find(ismember(zmos,seasons(k,:))==1);
                if ~isempty(seas_refs)
                    date_cat = ['    '];
                    for p = 1:length(seas_refs)
                        if p ~= 1
                            date_cat = [date_cat,', ',dateout(seas_refs(p),:)];
                        else
                            date_cat = [date_cat,dateout(seas_refs(p),:)];
                        end
                    end
                    disp(date_cat)
                end
                clear seas_refs;
            end

            %near-terminus thickness
            % disp(['    thickness (m): ',num2str(round(nanmean(MP(geo_ind(j)).B.Ho(1,k,:)),1))]);
            H(i,j,k) = nanmean(MP(geo_ind(j)).B.Ho(1,k,:));
            subplot(subH); %plot sensivity test results
            plot(MP(geo_ind(j)).Z.dist/10^3,MP(geo_ind(j)).Z.Hseas(k,:),'-','linewidth',2,'color',zcut_cmap(i,:)); hold on;
            if i == 3 %brute-force coding for the "best" zcutoff
                subplot(subHs);
                plot(MP(geo_ind(j)).Z.dist/10^3,MP(geo_ind(j)).Z.Hseas(k,:),'-','linewidth',2,'color',seas_cmap(k,:)); hold on;
            end

            %buttressing
            % disp(['    packing-based buttressing (N/m): ',num2str(round(nanmean(MP(geo_ind(j)).B.butt_Meng(1,k,:))./10^6,2))]);
            % disp(['    strainrate-based buttressing (N/m): ',num2str(round(nanmean(MP(geo_ind(j)).B.butt_Amundson(1,k,:))./10^6,2))]);
            bM(i,j,k) = nanmean(MP(geo_ind(j)).B.butt_Meng(1,k,:))./10^6;
            bA(i,j,k) = nanmean(MP(geo_ind(j)).B.butt_Amundson(1,k,:))./10^6;
            axes(ax1); %plot sensivity test results
            scatter(100*squeeze(MP(geo_ind(j)).B.packing(1,k,:)),squeeze(MP(geo_ind(j)).B.Ho(1,k,:)),15+15*squeeze(MP(geo_ind(j)).B.butt_Meng(1,k,:))./10^6,zcut_cmap(i,:),'o','LineWidth',1.5); hold on;
            axes(ax2);
            scatter(365*squeeze(MP(geo_ind(j)).B.dVdx(1,k,:)),squeeze(MP(geo_ind(j)).B.Ho(1,k,:)),15+15*squeeze(MP(geo_ind(j)).B.butt_Amundson(1,k,:))./10^6,zcut_cmap(i,:),'x','LineWidth',1.5); hold on;
            if i == 3 %brute-force coding for the "best" zcutoff
                axes(ax1s);
                scatter(100*squeeze(MP(geo_ind(j)).B.packing(1,k,:)),squeeze(MP(geo_ind(j)).B.Ho(1,k,:)),15+15*squeeze(MP(geo_ind(j)).B.butt_Meng(1,k,:))./10^6,seas_cmap(k,:),'o','LineWidth',1.5); hold on;
                axes(ax2s);
                scatter(365*squeeze(MP(geo_ind(j)).B.dVdx(1,k,:)),squeeze(MP(geo_ind(j)).B.Ho(1,k,:)),15+15*squeeze(MP(geo_ind(j)).B.butt_Amundson(1,k,:))./10^6,seas_cmap(k,:),'x','LineWidth',1.5); hold on;
            end
        end

        clear zdate datest dateout zyrs zmos;
    end

    clear MP;
end
%format the overlain sensitivity test subplots
ax2.XAxisLocation = 'top'; ax2.Color = 'none'; ax2.Box = 'on';
ax2.XLim = 365*[-0.0025,0.0225]; ax2.XTick = 365*[-0.0025:0.0025:0.0225];
ax2.XTickLabel = round(365*[-0.0025:0.0025:0.0225],2); Hylims = get(ax2,'ylim');
ax1.FontSize = 16; ax2.FontSize = 16;
ax2.XLabel.String = 'Strain rate (yr^{-1})'; ax2.XLabel.FontSize = 16;
ax1.XLabel.String = 'Packing density (%)'; ax1.XLabel.FontSize = 16;
ax1.YLabel.String = 'Thickness (m)'; ax1.YLabel.FontSize = 16;
%format the overlain seasonal subplots for the best zcutoff
ax2s.XAxisLocation = 'top'; ax2s.Color = 'none'; ax2s.Box = 'on';
ax2s.XLim = 365*[-0.0025,0.0225]; ax2s.XTick = 365*[-0.0025:0.0025:0.0225];
ax2s.XTickLabel = round(365*[-0.0025:0.0025:0.0225],2);
ax1s.FontSize = 16; ax1s.YLim = Hylims; ax2s.FontSize = 16; ax2s.YLim = Hylims;
ax2s.XLabel.String = 'Strain rate (yr^{-1})'; ax2s.XLabel.FontSize = 16;
ax1s.XLabel.String = 'Packing density (%)'; ax1s.XLabel.FontSize = 16;
ax1s.YLabel.String = 'Thickness (m)'; ax1s.YLabel.FontSize = 16;
%adjust the thickness profile plots
subplot(subH); grid on;
set(gca,'fontsize',16); 
ylabel('Thickness (m)','fontsize',16); xlabel('Distance from terminus (km)','fontsize',16); 
pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2)+0.05 pos(3) pos(4)]);
axpos = get(ax1,'position');
set(subH,'position',[0.15 axpos(2) axpos(3) axpos(4)]);
ylims = get(subH,'ylim');
axpos = get(ax1s,'position');
subplot(subHs); grid on;
set(subHs,'position',[0.15 axpos(2) axpos(3) axpos(4)]);
set(subHs,'ylim',ylims,'fontsize',16);
ylabel('Thickness (m)','fontsize',16); xlabel('Distance from terminus (km)','fontsize',16); 

%add varying legends & save
%colors legend for zcutoff sensitivity tests
subplot(subH); ylims = get(gca,'ylim');
rectangle('Position',[29,(0.98-0.04*6)*max(ylims),10.5,(0.04*6)*max(ylims)],'FaceColor','w');
for i = 1:5
    prefix = split(MPfiles(i).name,'zcutoff'); temp_name = char(prefix(1));
    plot([29.5 31],[(0.98-0.04*i)*max(ylims), (0.98-0.04*i)*max(ylims)],'-','linewidth',2,'color',zcut_cmap(i,:)); hold on;
    text(31.5,(0.98-0.04*i)*max(ylims),[temp_name(end-2:end-1),' z-cutoff'],'fontsize',14);
end
%colors legend for seasons
subplot(subHs); ylims = get(gca,'ylim');
rectangle('Position',[32,(0.98-0.04*5)*max(ylims),7.5,(0.04*5)*max(ylims)],'FaceColor','w');
for i = 1:4
    plot([32.5 34],[(0.98-0.04*i)*max(ylims), (0.98-0.04*i)*max(ylims)],'-','linewidth',2,'color',seas_cmap(i,:)); hold on;
    text(34.5,(0.98-0.04*i)*max(ylims),char(season_names(i)),'fontsize',14);
end
%shape legend for strain rates and packing density
axes(ax1); ylims = get(gca,'ylim'); grid on;
rectangle('Position',[22+2,0.88*max(ylims),30,0.1*max(ylims)],'FaceColor','w');
scatter(22+4,0.955*max(ylims),15+15,'k','o','LineWidth',1.5); hold on; text(22+6,0.955*max(ylims),'packing density','fontsize',14);
scatter(22+4,0.915*max(ylims),15+15,'k','x','LineWidth',1.5); hold on; text(22+6,0.915*max(ylims),'strain rate','fontsize',14);
rectangle('Position',[64+2,0.77*max(ylims),27,3*0.07*max(ylims)+0.02],'FaceColor','w');
for i = 1:3
    scatter(65+4,(0.94-(i-1)*0.07)*max(ylims),15+15*0.01*(10^i),'k','o','LineWidth',1.5); hold on;
    scatter(70+4,(0.94-(i-1)*0.07)*max(ylims),15+15*0.01*(10^i),'k','x','LineWidth',1.5); hold on;
    text(72+6,(0.95-(i-1)*0.07)*max(ylims),['10^',num2str(4+i),' N/m'],'fontsize',14);
end
axes(ax2);
%shape legend for buttressing based on symbol size
axes(ax1s); grid on;
rectangle('Position',[22+2,0.88*max(ylims),30,0.1*max(ylims)],'FaceColor','w');
scatter(22+4,0.955*max(ylims),15+15,'k','o','LineWidth',1.5); hold on; text(22+6,0.955*max(ylims),'packing density','fontsize',14);
scatter(22+4,0.915*max(ylims),15+15,'k','x','LineWidth',1.5); hold on; text(22+6,0.915*max(ylims),'strain rate','fontsize',14);
rectangle('Position',[64+2,0.77*max(ylims),27,3*0.07*max(ylims)+0.02],'FaceColor','w');
for i = 1:3
    scatter(65+4,(0.94-(i-1)*0.07)*max(ylims),15+15*0.01*(10^i),'k','o','LineWidth',1.5); hold on;
    scatter(70+4,(0.94-(i-1)*0.07)*max(ylims),15+15*0.01*(10^i),'k','x','LineWidth',1.5); hold on;
    text(72+6,(0.95-(i-1)*0.07)*max(ylims),['10^',num2str(4+i),' N/m'],'fontsize',14);
end
axes(ax2s);
saveas(zcutfig,[root_dir,'GrIS-melange-characteristics_sensitivity-plots.png'],'png'); %save the plots
exportgraphics(zcutfig,[root_dir,'GrIS-melange-characteristics_sensitivity-plots.tif'],Resolution=600);

%display the statistics for all seasons at each site
for j = 1:length(geo_ind)
    %site abbreviation and name
    disp(char(geo_names(j)));

    for k = 1:4
        if ~isnan(H(3,j,k))
            disp([' ',char(season_names(k))]);
            disp(['    thickness (m) = ',num2str(round(H(3,j,k),1)),' (',num2str(round(min(H(:,j,k)),1)),'-',num2str(round(max(H(:,j,k)),1)),')']);
            disp(['    packing-based buttressing (10^6 N/m) = ',num2str(round(bM(3,j,k),2)),' (',num2str(round(min(bM(:,j,k)),2)),'-',num2str(round(max(bM(:,j,k)),2)),')']);
            disp(['    strainrate-based buttressing (10^6 N/m) = ',num2str(round(bA(3,j,k),2)),' (',num2str(round(min(bA(:,j,k)),2)),'-',num2str(round(max(bA(:,j,k)),2)),')']);
        end
    end

end


%display relative differences in thickness & buttressing for the various
%zcutoff values
disp(['Fractional variability in thickness: ']);
for k = 1:4
    disp([char(season_names(k))]);
    round(nanmean(H(:,:,k),2)./nanmean(H(3,:,k),2),2)
end
disp(['Fractional variability in packing-based buttressing: ']);
for k = 1:4
    disp([char(season_names(k))]);
    round(nanmean(bM(:,:,k),2)./nanmean(bM(3,:,k),2),2)
end
disp(['Fractional variability in strain-based buttressing: ']);
for k = 1:4
    disp([char(season_names(k))]);
    round(nanmean(bA(:,:,k),2)./nanmean(bA(3,:,k),2),2)
end


