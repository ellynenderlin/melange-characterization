%% Create gifs of satellite images
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/');
% download S2 images with preautorift then make gifs
%identfy the files (run this section for each study site, specifying
%site-specific parameters below before each rerun)

%site-specific info
% root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';
root_dir = '/Users/ellynenderlin/Research/NASA_CryoIdaho/glaciers/';
site_abbrev = 'Wolverine'; site_name = 'Wolverine';
% im_dir = [root_dir,site_abbrev,'/images/S2/'];
im_dir = [root_dir,site_abbrev,'/imagery/S2/'];
disp(['Creating S2 image gif for ',site_abbrev]);

%adjust the time separation between images as needed (will crash if memory is exceeded!)
if strmatch('HLG',site_abbrev)
    ddays = 10; %at least 10 days between images
else
    ddays = 7; %at least a week between images
end

%load the data
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


    if contains(ims(k).name,'_2019') || str2num(im_dates(k,1:4)) >= 2020
        im_refs = [im_refs; k]; 
    end
end

%sort the images by date
for k = 1:length(im_refs)
    decidate(k,1) = convert_to_decimaldate(im_dates(im_refs(k),:)); 
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
xlims = get(gca,'xtick'); ylims = get(gca,'ytick');
xticks = get(gca,'xtick'); yticks = get(gca,'ytick');
set(gca,'xticklabels',xticks/1000,'yticklabels',yticks/1000,'fontsize',16);
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
title(site_name);
clear im;

%loop through the images and create a gif
nimages = 1;
for j = 1:length(date_refs)
    %decide whether to plot the image
    if j == 1
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
        
        frame = getframe(map_fig);
        gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    else
        if sorted_dates(j)-last_date >= ddays/365
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
            
            frame = getframe(map_fig);
            gif_im{nimages} = frame2im(frame); nimages = nimages+1;
        end
    end
    
    %     cla;
    %     title('');
    %     frame = getframe(map_fig);
    %     gif_im{nimages} = frame2im(frame); nimages = nimages+1;
    clear im;
end
close;

filename = [im_dir,site_abbrev,'-Sentinel2-images.gif']; % Specify the output file name
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