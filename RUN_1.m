% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% INPUTS
% Adminstrative Boundaries and Trajectory Parameters
%iso_3166_2 = {'US-CA','US-FL','US-KS','US-MA','US-MS','US-NC','US-ND','US-NH','US-NY','US-NV','US-OK','US-PR','US-RI','US-TN','US-TX','US-VA'}; % dev cases
iso_3166_2 = {'US-KS','US-MA','US-MS','US-NC','US-ND','US-NV','US-TN','US-VA'};
%iso_3166_2 = {'US-RI'}; % dev cases

file_airspace = [getenv('AEM_DIR_CORE') filesep 'output' filesep 'airspace-B-C-D-03-Aug-2020.mat'];

usecase = 'shield';

% Pairs parameters
% 30*(150 / 3600) = 1.25...closing speed of 150 knots for 30 seconds
% 2000 ft = well clear
anchorRange_nm = 1.00; % Distance pairs need to be away from anchor point 

%% Create status table to record performance
if ~iscolumn(iso_3166_2)
    iso_3166_2 = iso_3166_2';
end
status = table(iso_3166_2,zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),'VariableNames',{'iso_3166_2','runtime_s','n_points','n_files','n_clusters'});

%% Iterate through adminstrative boundaries
for i=1:1:numel(iso_3166_2)
    tic
    % Load Natural Earth Adminstrative Boundaries
    ne_admin = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'NE-Adminstrative' filesep 'ne_10m_admin_1_states_provinces.shp']);
    
    % Create bounding box
    BoundingBox_wgs84 = ne_admin(strcmp({ne_admin.iso_3166_2},iso_3166_2{i})).BoundingBox;
    
    % Load airspace
    load(file_airspace,'airspace');
    
    % Filter to include airspace only in bounding box
    [~, ~, inAirspace] = filterboundingbox(airspace.LAT_deg,airspace.LON_deg,BoundingBox_wgs84);
    airspace = airspace(inAirspace,:);
    
    neLat_deg = ne_admin(strcmp({ne_admin.iso_3166_2},iso_3166_2{i})).Y;
    neLon_deg = ne_admin(strcmp({ne_admin.iso_3166_2},iso_3166_2{i})).X;
    
    % Create grid
    [gridLat_deg, gridLon_deg] = meshgrid(min(neLat_deg):nm2deg(anchorRange_nm):max(neLat_deg),min(neLon_deg):nm2deg(anchorRange_nm):max(neLon_deg));
    
    % Filter grid for those only in the location
    isGrid = InPolygon(gridLat_deg,gridLon_deg,neLat_deg,neLon_deg);
    
    gridLat_deg = gridLat_deg(isGrid);
    gridLon_deg = gridLon_deg(isGrid);
    
    % Uncomment to plot
    % figure(i); set(gcf,'name',iso_3166_2{i}); worldmap(BoundingBox_wgs84(:,2),BoundingBox_wgs84(:,1)); geoshow(gridLat_deg,gridLon_deg,'DisplayType','point'); grid on;
    
    % Input directories for ownship files
    % The following code assumes you're using em-model-geospatial.
    % If not, replace this section with code to identify the filenames of your ownship trajectory files
    listing = dir([getenv('AEM_DIR_GEOSPATIAL') filesep 'output' filesep 'trajectories' filesep iso_3166_2{i} '*' filesep '*']);
    listing(ismember( {listing.name}, {'.', '..'})) = [];  %remove . and ..
    isShield = contains({listing.folder},'shield');
    inDir = cellfun(@(f,n)([f filesep n]),{listing.folder},{listing.name},'UniformOutput',false)';
    switch usecase
        case 'all'
            % No filtering needed
        case 'all-noshield'
            inDir(isShield) = [];
        case 'unconv-noshield'
            % Pairs to be used with the unconventional bayes model
            %(i.e. paragliders, gliders,)
            isCostal = contains({listing.folder},{'landuse_beach','landuse_cliff','landuse_volcano','gshhg_gshhs_resf_lvl1'});
            inDir(~isCostal) = [];
        case 'shield'
            inDir(~isShield) = [];
        case 'longlinearinfrastructure'
            islli = contains({listing.folder},{'pipeline','roads','waterway','railway','electrictransmission'});
            inDir(~islli) = [];
        case 'railway'
            irrail = contains({listing.folder},'railway');
            inDir(~irrail) = [];
        case 'agriculture'
            isAg = contains({listing.folder},{'farm','orchard','vineyard'});
            inDir(~isAg) = [];
        case 'heliport5-all-noshield'
            inDir(isShield) = [];
            % Just filtering for grid coordinates
            % Load airport data and filter to heliports
            S_airports = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'FAA-Airports' filesep 'Airports.shp'],'BoundingBox',BoundingBox_wgs84,'UseGeoCoords',true);
            S_hp = S_airports(strcmpi({S_airports.TYPE_CODE},'HP'));
            
            % Determine which points are within 5 nm of airport
            inHP = sparse(numel(gridLat_deg),size(S_hp,1));
            for j=1:1:size(S_hp,1)
                [latc,lonc] = scircle1(S_hp(j).Lat,S_hp(j).Lon,5,[],wgs84Ellipsoid('nm'));
                inHP(:,j) = InPolygon(gridLat_deg,gridLon_deg,latc,lonc);
            end
            % Aggregate across rows
            inHP = any(inHP,2);
            
            % Filter
            gridLat_deg = gridLat_deg(find(inHP));
            gridLon_deg = gridLon_deg(find(inHP));
            
            % Plot
            figure(i); set(gcf,'name',iso_3166_2{i});
            gx = geoaxes('Basemap','streets-dark',...
                'FontSize',16,...
                'FontWeight','bold',...
                'Position',[0 0 1 1]);
            bmap = geobasemap(gx); hold on;
            geoscatter(gridLat_deg,gridLon_deg,'o','MarkerEdgeColor',[86 180 233]/255);
            geoscatter([S_hp.Lat],[S_hp.Lon],'s','filled','MarkerFaceColor',[240 228 66]/255);
            geolimits( [BoundingBox_wgs84(1,2) - nm2deg(6),  BoundingBox_wgs84(2,2) + nm2deg(6)], [BoundingBox_wgs84(1,1) - nm2deg(6),  BoundingBox_wgs84(2,1) + nm2deg(6)]);    hold off;
            %legend('Potential Anchor Points','Heliports');
            set(gcf,'Units','inches','Position',[1 1 11.94 5.28]); % Adjust map size
            print(['output' filesep 'plot_grid_' usecase '_' iso_3166_2{i} '.png'],'-dpng','-r300');
        otherwise
            error('usecase:unknown','Unknown use case of %s\n',usecase);
    end
    
    if isempty(inDir)
        fprintf('No input directories for usecase = %s, iso_3166_2 = %s...skipping via CONTINUE\n',usecase,iso_3166_2{i});
        continue;
    end
    
    % Execute
    [anchors, inFiles, numClust] = findPairs_1(inDir,sprintf('%s_%s',usecase,iso_3166_2{i}),gridLat_deg,gridLon_deg,airspace,'anchorRange_nm',anchorRange_nm);
    
    % Record stats
    status.runtime_s(i) = round(toc);
    status.n_points(i) = numel(gridLat_deg);
    status.n_files(i) = numel(inFiles);
    status.n_clusters(i) = numClust;
    
    % These are saved and we only outputted for performance, don't need
    % them hanging around for the next loop
    clear anchors inFiles
end

%% Save Display status to screen
save(['output' filesep '1_performanace_' usecase '_' 'anchorRange' num2str(anchorRange_nm) '_' date '.mat'],'status','inDir','file_airspace','anchorRange_nm','usecase');
disp('Done!');
