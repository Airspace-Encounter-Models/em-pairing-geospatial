% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [anchors, listing, numClust] = findPairs_1(inDir,outHash,gridLat_deg,gridLon_deg,airspace,varargin)

%% Things that should be input parser

%% Input parser
p = inputParser;

% Required
addRequired(p,'inDir',@iscell);
addRequired(p,'outHash');
addRequired(p,'gridLat_deg',@isnumeric);
addRequired(p,'gridLon_deg',@isnumeric);
addRequired(p,'airspace');

% Optional
addOptional(p,'anchorRange_nm',1.25,@isnumeric);
addOptional(p,'classInclude',{'B','C','D','O'},@iscell);

addOptional(p,'spheroid',wgs84Ellipsoid('nm'));
addOptional(p,'isSave',true,@islogical);

% Parse
parse(p,inDir,outHash,gridLat_deg,gridLon_deg,airspace,varargin{:});
anchorRange_nm = p.Results.anchorRange_nm;

%% Create outputs 
% table
anchors = table(gridLat_deg,gridLon_deg,repmat('O',size(gridLon_deg)),cell(size(gridLon_deg)),zeros(size(gridLon_deg)),'VariableNames',{'lat_deg','lon_deg','class','files','num_geospatial'});

% filename
outName = ['output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat',outHash,anchorRange_nm)];

% Preallocate
numClust = 0;

%% Parse filenames and metadata
% Get list of filenames
listing = dir([inDir{1} filesep '*.csv']);

% Concat listings if there is more than one input directory
if numel(inDir) > 1
    for i=2:1:numel(inDir)
        listing = [listing; dir([inDir{i} filesep '*.csv'])];
    end
end

if isempty(listing); warning('RUN_1:listing','No files found inDir, exiting via RETURN \n'); return; end

% Set consistent random seed
rng(size(listing,1));

% Filter trajectories based on
isAir = false(size(listing,1),1);
minLon = zeros(size(listing,1),1);
minLat = zeros(size(listing,1),1);
maxLon = zeros(size(listing,1),1);
maxLat = zeros(size(listing,1),1);

parfor i=1:1:size(listing,1)
    C = strsplit(listing(i).name,'_');
    lclass = contains(C,'class','IgnoreCase',true);
    if any(lclass)
        isAir(i) = contains(C{lclass}(6:end),p.Results.classInclude);
    else
        % Not every model has airspace class, like HAA
        isAir(i) = true;
    end
    
    % Parse bounding box if airspace is accepted
    if isAir(i)
        llonmin = contains(C,'lonmin','IgnoreCase',true);
        if any(llonmin); minLon(i) = str2double(C{llonmin}(7:end)); end
        
        llonmax = contains(C,'lonmax','IgnoreCase',true);
        if any(llonmax); maxLon(i) = str2double(C{llonmax}(7:end)); end
        
        llatmin = contains(C,'latmin','IgnoreCase',true);
        if any(llatmin); minLat(i) = str2double(C{llatmin}(7:end)); end
        
        llatmax = contains(C,'latmax','IgnoreCase',true);
        if any(llatmax); maxLat(i) = str2double(C{llatmax}(7:end)); end;
    end
end
listing(~isAir) = [];
minLon(~isAir) = [];
minLat(~isAir) = [];
maxLon(~isAir) = [];
maxLat(~isAir) = [];

%% Create clusters
area_deg = abs(max(maxLon) - min(minLon)) * abs(max(maxLat) - min(minLat));
area_nm = deg2nm(area_deg);
numClust = ceil(area_nm / 100);
[idx,~,~] = kmeans([minLat,minLon,maxLat,maxLon],numClust);
fprintf('Iterating over %i kmeans clusters\n',numClust);

%% Iterate over clusters
for k=1:1:numClust
    % Filter listing to those in cluster
    isk = idx == k;
    listk = listing(isk);
    fprintf('Cluster %i: %i features\n', k, sum(isk));
    
    % Filter potential anchor points to bounding box of cluster
    buff_deg = anchorRange_nm;
    bbox = [min(minLon(isk))-buff_deg, min(minLat(isk))-buff_deg; max(maxLon(isk))+buff_deg, max(maxLat(isk))+buff_deg];
    
    % Filter candidate anchor points
    [candLat, candLon, isCluster] = filterboundingbox(gridLat_deg,gridLon_deg,bbox);
    
    % Display to screen and CONTINUE if needed
    if numel(candLat) > 0
        fprintf('Cluster %i: %i candidate anchor points\n',k, numel(candLat));
    else
        % This is possible because some of the geospatial representative
        % trajectories are filtered based on bounding box, not polygons.
        % For example, the US-CA bounding box will include multiple states,
        % whereas (gridLat_deg, gridLon_deg) could be based on the US-CA polygon
        fprintf('Cluster %i: %i candidate anchor points...calling CONTINUE\n',k, numel(candLat));
        continue;
    end
    
    % Create small circle centered on anchor point with radius of anchorRange_nm
    [latc,lonc] = scircle1(candLat,candLon,repmat(anchorRange_nm,numel(candLat),1),[],p.Results.spheroid);
    
    % Filter airspace based on bounding box
    [~, ~, inAir] = filterboundingbox(airspace.LAT_deg,airspace.LON_deg,bbox);
    
    % Filter airspace based on altitude surface
    isAir = inAir & cellfun(@min,airspace.LOWALT_ft_agl) <= 0;
    
    % Airspace
    airB = airspace(airspace.CLASS == 'B' & isAir,:);
    airC = airspace(airspace.CLASS == 'C' & isAir,:);
    airD = airspace(airspace.CLASS == 'D' & isAir,:);
    airE = airspace(airspace.CLASS == 'E' & isAir,:);
    
    % Identify airspace class
    % This is slow
    [isE,~] = identifyairspace(airE, candLat, candLon);
    [isD,~] = identifyairspace(airD, candLat, candLon);
    [isC,~] = identifyairspace(airC, candLat, candLon);
    [isB,~] = identifyairspace(airB, candLat, candLon);
    
    candAir = anchors.class(isCluster);
    candAir(isE) = 'E';
    candAir(isD) = 'D';
    candAir(isC) = 'C';
    candAir(isB) = 'B';
        
    % Get coordinates for each trajectory
    % Preallocate
    time_s = cell(size(listk));
    lat_deg = cell(size(listk));
    lon_deg = cell(size(listk));
    bboxlat = cell(size(listk));
    bboxlon = cell(size(listk));
    
    % Iterate over files
    parfor i=1:1:numel(listk)
        % Load trajectory
        [time_s{i}, lat_deg{i}, lon_deg{i}, ~] = fastloadtraj([listk(i).folder filesep listk(i).name])
        
        % Calculate bounding box
        lonmin = min(lon_deg{i});
        lonmax = max(lon_deg{i});
        latmin = min(lat_deg{i});
        latmax = max(lat_deg{i});
        
        bboxlon{i} = [lonmin lonmax lonmax lonmin lonmin];
        bboxlat{i} = [latmin latmin latmax latmax latmin];
    end
    
    % For candidate pairs, determine if they are in the bounding boxes
    l = cellfun(@(x,y)(InPolygon(candLon,candLat,x,y)),bboxlon,bboxlat,'uni',false);
    
    % For each bounding box, determine which of s:e pairs
    isPair = cellfun(@find,l,'uni',false);
    
    % Preallocate variable for parfor because inBox can't be used
    inBox = cell(numel(candLon),1);
    
    % Logical indx of trajectories that have no potential overlap with
    % candidate anchor points, we do this to speed up iterating through
    % each pair
    noPair = cellfun(@isempty,isPair);
    
    % Iterate through each pair
    parfor i=1:1:numel(inBox)
        % Preallocate
        inBox{i} = false(size(isPair));
        
        % Find if ith candidate for each trajectory
        inBox{i}(~noPair) = cellfun(@(x)(~isempty(find(x==i))),isPair(~noPair),'uni',true);
        %cellfun(@(x)(ismember(i,x)),isPair);
    end
    
    % Iterate through candidates
    candFiles = anchors.files(isCluster);
    for i=1:1:numel(candLat)
        
        % Get index for lat_deg, lon_deg for trajectories whose bounding box
        % includes the anchor point i
        idxTraj = find(inBox{i} == true);
        
        if any(idxTraj)
            % Calculate distance from anchor points
            isClose = cellfun(@(x,y)(any(InPolygon(x,y,lonc(:,i),latc(:,i)))),lon_deg(idxTraj),lat_deg(idxTraj),'uni',true);
            
            % Filter index to those only close to the anchor point
            idxTraj = idxTraj(isClose);
            
            % Create full filenames for trajectories which are close to anchor point
            files = compose('%s/%s',string({listk(idxTraj).folder}'), string({listk(idxTraj).name}'));
            
            % Assign
            candFiles{i} = [candFiles{i};files];
        end
        %
        %         [lat,lon] = polyjoin(lat_deg(idxTraj),lon_deg(idxTraj));
        %         geoshow(lat,lon); hold on
        %         geoshow(latc(:,i),lonc(:,i));
    end
    
    % Assign to output
    anchors.class(isCluster) = candAir;
    anchors.files(isCluster) = candFiles; 
end % End k loop

%% Clean and Finish Up
% Remove duplicate files
anchors.files = cellfun(@unique,anchors.files,'UniformOutput',false);

% Convert to string
anchors.files = cellfun(@string,anchors.files,'uni',false);

% Calculate number of trajectories for each anchor point
anchors.num_geospatial = cellfun(@numel,anchors.files);

%% Save to file
if p.Results.isSave
    save(outName,'anchors','anchorRange_nm');
end
