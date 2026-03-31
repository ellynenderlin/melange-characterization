function [berg_dates,Havg,packing,inland_idx,seaward_idx,seaward_ext,H_seas,pack_seas] = extract_seasonal_binned_DEM_profiles(D_dir,site_abbrev,term_trace,term_x,term_y,center_x,center_y,WHratio,zcutoff,sampling,years,seasons)
%%% Extract seasonal elevation-based profiles from binned iceberg size distributions

%INPUTS:
%D_dir = binned size distribution CSV directory
%site_abbrev = name of site that is at the start of the CSV files
%term_trace = flag indicating if a terminus was visible in a DEM (1) or interpolated from an image with a similar date (0)
%inland_idx = index or the point along the profile where the melange begins
%term_x,term_y = x and y coordinates for the terminus intersection with the profile from satellite images
%center_x,center_y = x and y coordinates for the profile from which the elevation bins were defined
%WHratio = iceberg width-to-thickness ratio for converting elevations to size distributions (needed to back-convert to elevation distributions)
%zcutoff = elevation threshold (z < zcutoff are bergy bits not included in calculations)
%sampling = flag to determine whether to use the same profile point as the reference for all data ('fixed') or based on the terminus position ('date')
%years = years that you want in your speed dataset
%seasons = months corresponding to each season (recommend: seasons = [12,1,2;3,4,5;6,7,8;9,10,11]; season_names = {'DJF','MAM','JJA','SON'};)

%OUTPUTS:
%berg_dates = date strings for elevation data (yyyymmdd)
%Havg = average thickness profile for each date (m)
%packing = average packing density profile for each date (m)
%inland_idx = point reference for inland data extent 
%seaward_idx = point reference for TYPICAL seaward data extent (used for plotting size distributions)
%seaward_ext = point reference for seaward data extent
%H_seas = average thickness profile for each season (m)
%pack_seas = average packing density profile for each season (m)

%check that these constants match what was used to convert elevations to
%thicknesses in "extract_automated_iceberg_DEM_distributions.m"
rho_i = 900; rho_sw = 1026; %density of ice and sea water in kg/m^3 (constant)

%loop through subset size distributions & find the first and last
%non-NaN columns to identify melange extent
Dsubs = dir([D_dir,site_abbrev,'*-iceberg-distribution-subsets.csv']);
size_classes = [];
for p = 1:length(Dsubs)
    D = readtable([D_dir,Dsubs(p).name],"VariableNamingRule","preserve");
    name_split = split(Dsubs(p).name,'-',2);
    berg_dates(p) = datetime(name_split{2},'Inputformat','yyyyMMdd');
    zdate(p) = convert_to_decimaldate(char(name_split{2})); clear name_split;
    zyrs(p) = year(berg_dates(p)); zmos(p) = month(berg_dates(p));

    %identify observational limits along the centerline
    berg_nos = table2array(D(:,3:end)); berg_nos(berg_nos==0) = NaN;
    berg_sizes = ((rho_sw-rho_i)/rho_sw)*((2/WHratio)*sqrt(table2array(D(:,1))/pi())); dsize = round(nanmean(diff(berg_sizes)));
    berg_maxz = berg_sizes+0.5*dsize; size_ind = find(berg_maxz <= zcutoff,1,'last');
    size_classes(p,:) = sum(~isnan(berg_nos),1);
    seaward_ext(p) = find(size_classes(p,:)>0,1,'first');
    inland_ext(p) = find(size_classes(p,:)>0,1,'last')+1;

    %create average thickness profile
    Havg(p,:) = sum((2/WHratio)*sqrt(table2array(D(size_ind+1:end,1))./pi()).*table2array(D(size_ind+1:end,3:end)),1)./sum(table2array(D(size_ind+1:end,3:end)),1);
    packing(p,:) = sum(table2array(D(size_ind+1:end,1)).*table2array(D(size_ind+1:end,3:end)),1)./sum(table2array(D(:,1)).*table2array(D(:,3:end)),1);
    clear size_ind;
end

%assign seaward & inland sampling limits based on the selected sampling
%strategy (sampling = ['fixed','dated']
inland_meanidx = round(mean(inland_ext(term_trace==1)));
seaward_meanidx = round(mean(seaward_ext));
if contains(sampling,'fix')
    inland_idx = inland_meanidx.*ones(size(seaward_ext));
    seaward_idx = seaward_meanidx.*ones(size(seaward_ext));
elseif contains(sampling,'date')
    %if the terminus was visible in the DEM, use it as your reference
    inland_idx = inland_ext;

    %if the terminus was NOT visible in the DEM, use a combination of
    %the DEM terminus timeseries and TermPicks timeseries to
    %approximate the terminus location (gaps were filled for termini
    %with observations within 60 days in the previous section)
    for p = 1:length(term_trace)
        if term_trace(p) == 0
            if ~isnan(term_x(p))
                dists = sqrt((term_x(p)-center_x).^2 + (term_y(p)-center_y).^2);
                [~,ia] = sort(dists);
                inland_idx(p) = max(ia(1:2)); term_trace(p) = 1; %use terminus traces from imagery
                clear dists ia;
            else
                inland_idx(p) = NaN;
                disp(['Need terminus data for ',site_abbrev,' ',char(berg_dates(p))])
            end
        end
    end
    % seaward_idx = inland_idx-(inland_meanidx-seaward_meanidx);
    seaward_idx = inland_idx-floor(nanmean(inland_idx-seaward_ext));
else
    error('Profile sampling strategy is not properly defined: sampling must be ''fixed'' or ''dated''')
end

%create an annual average seasonal thickness profile
% Zfilt = MP(j).Z.transectZavg; Zfilt(MP(j).Z.transectZavg==0) = NaN; Zfilt(:,term_trace==0) = NaN;
Havg(term_trace==0,:) = NaN; H_seas = NaN(max(inland_idx)-1,4,length(years)); pack_seas = NaN(max(inland_idx)-1,4,length(years));
for p = 1:length(years)
    yr_idx = find(zyrs == years(p));
    if ~isempty(yr_idx)
        mos_yr = zmos(yr_idx);

        %create a matrix of elevations relative to the terminus position
        % z_profiles = NaN(max(inland_idx),length(yr_idx));
        H_profiles = NaN(length(yr_idx),max(inland_idx)-1);
        pack_profiles = NaN(length(yr_idx),max(inland_idx)-1);
        for k = 1:length(yr_idx)
            if ~isnan(inland_idx(yr_idx(k)))
                % z_temp = flipud(Zfilt(1:inland_idx(yr_idx(k))-1,yr_idx(k)));
                % z_profiles(1:length(z_temp),k) = z_temp;
                H_temp = fliplr(Havg(yr_idx(k),1:inland_idx(yr_idx(k))-1));
                H_profiles(k,1:length(H_temp)) = H_temp;
                pack_temp = fliplr(packing(yr_idx(k),1:inland_idx(yr_idx(k))-1));
                pack_profiles(k,1:length(pack_temp)) = pack_temp;
                clear H_temp pack_temp; %clear z_temp;
            end
        end

        %calculate seasonal averages
        for k = 1:4
            % z_seas(:,k,p) = nanmean(z_profiles(:,ismember(mos_yr,seasons(k,:))==1),2);
            H_seas(:,k,p) = nanmean(H_profiles(ismember(mos_yr,seasons(k,:))==1,:),1)';
            pack_seas(:,k,p) = nanmean(pack_profiles(ismember(mos_yr,seasons(k,:))==1,:),1)';
        end

        clear H_profiles pack_profiles mos_yr; %clear z_profiles;
    end
    clear yr_idx
end
% z_seas(z_seas==0) = NaN;
H_seas(H_seas==0) = NaN; pack_seas(pack_seas==0) = NaN;

end