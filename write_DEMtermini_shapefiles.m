%%% Create DEM-based terminus position shapefiles
clearvars;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/',...
    '/Users/ellynenderlin/Research/miscellaneous/general-code/inpoly2/');
addpath('/Users/ellynenderlin/Research/miscellaneous/melange-characterization-code/');

%specify root path
root_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/melange/';

%load the data
cd(root_dir);
sites = dir; sitenames = [];
for i = 1:length(sites)
    if ~contains(sites(i).name,'.') && length(sites(i).name) == 3
        sitenames = [sitenames; sites(i).name];
    end
end
load([root_dir,'GrIS-melange-characteristics.mat']);

%iterate
for j = 1:length(MP)
    disp(sitenames(j,:)); output_dir = [root_dir,sitenames(j,:),'/termini/']; site_abbrev = MP(j).name;

    %load the melange masks
    cd([root_dir,sitenames(j,:)]);
    load([MP(j).name,'-melange-masks.mat']); %load the melange mask file

    %isolate terminus positions from the dated melange masks
    % T = geoshape();
    for p = 1:length(melmask.dated)
        [~,on] = inpolygon(melmask.dated(p).x,melmask.dated(p).y,melmask.uncropped.x,melmask.uncropped.y);
        melmask.dated(p).x(on) = []; melmask.dated(p).y(on) = [];
        % T = append(T,[melmask.dated(p).y; NaN],[melmask.dated(p).x; NaN],'Date',melmask.dated(p).datestring);
        Ttable(p).Date = {melmask.dated(p).datestring};
        Ttable(p).X = [melmask.dated(p).x; NaN]'; Ttable(p).Y = [melmask.dated(p).y; NaN]';
    end
    %write the data as a geotable with the correct projection
    crs = projcrs(3413); 
    Tgeotable = struct2geotable(Ttable,CoordinateReferenceSystem=crs);
    shapewrite(Tgeotable,[root_dir,sitenames(j,:),'/termini/',MP(j).name,'-ArcticDEMtermini_2011-2023.shp']);
    
    %override the default "mappointshape" written for the geotable so each
    %terminus delineation is a polyline
    S = shaperead([root_dir,sitenames(j,:),'/termini/',MP(j).name,'-ArcticDEMtermini_2011-2023.shp']);
    for p = 1:length(S)
        S(p).Geometry = 'PolyLine';
    end
    shapewrite(S,[root_dir,sitenames(j,:),'/termini/',MP(j).name,'-ArcticDEMtermini_2011-2023.shp']);
    disp(['DEM terminus delineation shapefile written for ',MP(j).name]);

    clear melmask on T*table S;
end

