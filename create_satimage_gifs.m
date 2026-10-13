%% Create gifs of satellite images
% download S2 images with preautorift then make gifs
clearvars; close all; warning off;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/');

%% Define paths & gif details

%site-specific info
% root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-melt-data/images/';
site_abbrev = ['ILS']; 
% im_dir = [root_dir,site_abbrev,'/images/S2/'];
im_dir = [root_dir,site_abbrev,'/'];
disp(['Creating S2 image gif for ',site_abbrev]);

% specify the purpose of the gifs
gif_purpose = 'icebergs'; 
%options: 'icebergs' will add symbols for iceberg locations, 
% 'melange' will grab melange images from a different folder but not show icebergs, 
% 'general' will just show the images with no overlays
if contains(gif_purpose,'icebergs') || contains(gif_purpose,'melange')
    melange_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
    melange_info = dir(melange_dir); melange_dirs = [];
    for k = 1:length(melange_info)
        if isfolder([melange_dir,melange_info(k).name]) == 1 && length(melange_info(k).name) == 3
            melange_dirs = [melange_dirs; melange_info(k).name];
        end
    end

    %point to images compiled for the melange for im_dir
    if sum(contains(string(melange_dirs),site_abbrev)) >= 1
        im_dir = [melange_dir,site_abbrev,'/images/S2/'];
    end
end

%adjust the time separation between images as needed (will crash if memory is exceeded!)
if strmatch('HLG',site_abbrev)
    ddays = 10; %at least 10 days between images
else
    ddays = 7; %at least a week between images
end

%% Create the gif

% load the coordinates for icebergs if specified
if contains(gif_purpose,'icebergs')
    %read the DEM datepair dates and X,Y positions from the site CSV
    melt_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-melt-data/site-csvs/';
    melt = readtable([melt_dir,site_abbrev,'_basic_iceberg_meltinfo.csv']);
    for k = 1:size(melt,1)
        start_date = table2array(melt(k,2)); end_date = table2array(melt(k,3));
        start_datetime(k) = datetime(num2str(start_date),'InputFormat','yyyyMMdd');
        end_datetime(k) = datetime(num2str(end_date),'InputFormat','yyyyMMdd');
        PS_X(k) = melt.X_m_(k); PS_Y(k) = melt.Y_m_(k);
        [lon,lat] = ps2wgs(PS_X(k),PS_Y(k));
        [X(k),Y(k),~,~] = wgs2utm(lat,lon);

        clear start_date end_date lat lon;
    end
end

%load the image dates & sort
% load([root_dir,site_abbrev,'/',site_abbrev,'-melange-masks.mat']); %load the melange mask file
% ims = dir([im_dir,'S*B08_clipped.tif']); im_refs = []; im_dates = [];
ims = dir([im_dir,'S*_clipped.tif']); im_refs = []; im_dates = [];
for k = 1:length(ims)
    name_split = split(ims(k).name,'_',2);
    if length(char(name_split(2))) ==5
        im_dates = [im_dates; ims(k).name(11:18)];
    else
        im_dates = [im_dates; ims(k).name(10:17)];
    end

    %filter dates by year
    % if contains(ims(k).name,'_2019') || str2num(im_dates(k,1:4)) >= 2020
    if str2num(im_dates(k,1:4)) >= 2017
        im_refs = [im_refs; k]; 
    end
end
%sort the images by date
for k = 1:length(im_refs)
    decidate(k,1) = convert_to_decimaldate(im_dates(im_refs(k),:)); 
    im_datetime(k) = datetime(im_dates(im_refs(k),:),'InputFormat','yyyyMMdd');
end
[sorted_dates,date_refs] = sort(decidate);

%load the first image to serve as a reference 
% [I,R] = readgeoraster([im_dir,ref_image]);
[I,R] = readgeoraster([im_dir,ims(im_refs(date_refs(1))).name]);
im.x = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
im.y = linspace(R.YWorldLimits(2),R.YWorldLimits(1),R.RasterSize(1));
im.z = double(I);
clear I R;

%create the map template with the reference date shown to start
map_fig = figure; set(map_fig,'position',[850 50 800 600]);
imagesc(im.x,im.y,imadjust(im.z./max(max(im.z)))); axis xy equal; colormap gray; drawnow; hold on;
if contains(site_abbrev,'SEK')
    set(gca,'xlim',[537000 559000],'ylim',[7668000 7683000]); 
else
    set(gca,'xlim',[min(im.x), max(im.x)],'ylim',[min(im.y),max(im.y)]);
end
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
clear im;

%loop through the images and create a gif
nimages = 1;
for j = 1:length(date_refs)
    %decide whether to plot the image
    if j == 1 || sorted_dates(j)-last_date >= ddays/365
        %load the image for the specified date
        [I,R] = readgeoraster([im_dir,ims(im_refs(date_refs(j))).name]);
        im.x = linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
        im.y = linspace(R.YWorldLimits(2),R.YWorldLimits(1),R.RasterSize(1));
        im.z = double(I);
        clear I R;
        %plot
        imagesc(im.x,im.y,imadjust(im.z./max(max(im.z)))); axis xy equal; colormap gray; drawnow; hold on;
        if contains(site_abbrev,'SEK')
            set(gca,'xlim',[537000 559000],'ylim',[7668000 7683000]);
        else
            set(gca,'xlim',[min(xlims),max(xlims)],'ylim',[min(ylims),max(ylims)]);
        end
        set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
        xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
        title([im_dates(im_refs(date_refs(j)),1:4),'/',im_dates(im_refs(date_refs(j)),5:6),'/',im_dates(im_refs(date_refs(j)),7:8)]);
        drawnow;
        last_date = sorted_dates(j);

        % show the avg location of icebergs for observations between melange DEM date pairs
        if contains(gif_purpose,'icebergs')
            for k = 1:length(X)
                if isbetween(im_datetime(im_refs(date_refs(j))),start_datetime(k),end_datetime(k))
                    plot(X(k),Y(k),'*','MarkerEdgeColor',[1 0.5 0],'MarkerFaceColor',[1 0.5 0],'markersize',12,'linewidth',2); hold on;
                end
            end
        end

        frame = getframe(map_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    end
    
    %     cla;
    %     title('');
    %     frame = getframe(map_fig);
    %     gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    clear im;
end
close;

% define the gif output path and name
if contains(gif_purpose,'icebergs')
    filename = [root_dir,site_abbrev,'/',site_abbrev,'-icebergs_Sentinel2-images.gif']; % Specify the output file name
else
    filename = [root_dir,site_abbrev,'/',site_abbrev,'_Sentinel2-images.gif']; % Specify the output file name
end
% export the gif
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
close all; clear A frame gif_im map;
clear sort* im_* ims date_refs decidate;
disp('GIF written!');