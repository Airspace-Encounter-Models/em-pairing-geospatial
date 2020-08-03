% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [encCount, anchors, files1, files2, isEncounter] = createEncounters_2(inFile,encTime_s,thresHorz_ft,thresVert_ft,outDir,varargin)
% CREATEENCOUNTERS_2 OUTPUTS ENCOUNTERS IN CSIM WAYPOINT FORMAT

%% Input parser
p = inputParser;

% Required
addRequired(p,'inFile',@ischar); % Output of RUN_1
addRequired(p,'encTime_s',@isnumeric); % Duration of encounter
addRequired(p,'thresHorz_ft',@isnumeric);
addRequired(p,'thresVert_ft',@isnumeric);
addRequired(p,'outDir'); % Output directory

% Optional - Minimum encounter initial conditions
addOptional(p,'initHorz_ft',[6076 121522],@(x) isnumeric(x) && numel(x) == 2); %121522 ~ 20 nm
addOptional(p,'initVert_ft',[0 300],@(x) isnumeric(x) && numel(x) == 2);
addOptional(p,'spheroid',wgs84Ellipsoid('ft'));

% Optional - Anchor
addOptional(p,'classInclude',{'B','C','D','O'},@iscell);

% Optional - Encounter Variations
addOptional(p,'maxEncPerPair',10,@isnumeric); % Maximum waypoint pairs for each unique pair
addOptional(p,'maxCombo',1e7,@isnumeric); %  Maximum waypoint combinations to consider for each pair
addOptional(p,'timeStep_s',1,@(x) isnumeric(x) && (x >=1)); %  Output timestep, CSIM requires at least >= 1
addOptional(p,'anchorPercent',1.0,@isnumeric); %  Maximum waypoint combinations to consider for each pair

% Optional - AC permutations
addOptional(p,'acmodel1','geospatial',@(x) ischar(x) && any(strcmpi(x,{'geospatial'})));
addOptional(p,'acmodel2','geospatial',@(x) ischar(x) && any(strcmpi(x,{'bayes','geospatial'})));
addOptional(p,'inDir1',{},@iscell);
addOptional(p,'inDir2',{},@iscell);

% Optional - Airspeed Sampling (Geospatial only)
addOptional(p,'vrange1_kts',[5 87],@(x) isnumeric(x) && numel(x) == 2);
addOptional(p,'vrange2_kts',[5 87],@(x) isnumeric(x) && numel(x) == 2);

% Optional - Altitude Sampling (Geospatial only)
addOptional(p,'isSampleAlt1',false,@islogical);
addOptional(p,'isSampleAlt2',false,@islogical);
addOptional(p,'minAlt1_ft_agl',100,@isnumeric); % Minimum AGL altitude
addOptional(p,'minAlt2_ft_agl',100,@isnumeric); % Minimum AGL altitude
addOptional(p,'maxAlt1_ft_agl',400,@isnumeric); % Maximum AGL altitude
addOptional(p,'maxAlt2_ft_agl',400,@isnumeric); % Maximum AGL altitude

% Optional - Bayes Tracks
addOptional(p,'z_agl_tol_ft',200,@(x) isnumeric(x) && numel(x) == 1); % Used by placeTrack()
addOptional(p,'demBayes','strm3',@(x) isstr(x) && any(strcmpi(x,{'dted1','dted2','globe','gtopo30','srtm1','srtm3','srtm30'})));
addOptional(p,'demBayesBackup','globe',@(x) isstr(x) && any(strcmpi(x,{'dted1','dted2','globe','gtopo30','srtm1','srtm3','srtm30'})));
addOptional(p,'labelTime','time_s');
addOptional(p,'labelX','x_ft');
addOptional(p,'labelY','y_ft');
addOptional(p,'labelZ','z_ft');

% Optional - FAA DOF
addOptional(p,'dofFile',[getenv('AEM_DIR_CORE') filesep 'output' filesep 'dof-' '25-Feb-2020' '.mat'],@ischar); % Location of output of RUN_readfaadof from em-core
addOptional(p,'dofObs',{'grain elevator','lgthouse','pole','silo','tower','utility pole','windmill'},@iscell); % Obstacles to keep
addOptional(p,'dofMinHeight_ft',100,@isnumeric); % Minimum height of obstacles
addOptional(p,'dofMaxRange_ft',500,@isnumeric);
addOptional(p,'dofMaxVert_ft',50,@isnumeric);

% Optional - Output
addOptional(p,'isPlot',false,@islogical);
addOptional(p,'isZip',true,@islogical);

% Row offsets
% These are needed because the geospatial trajectory generator will add
% vertical take off and landing segements, which you may not want
addOptional(p,'rowstart1',2,@isnumeric); %
addOptional(p,'rowend1',1,@isnumeric); %
addOptional(p,'rowstart2',2,@isnumeric); %
addOptional(p,'rowend2',1,@isnumeric); %

% Parse
parse(p,inFile,encTime_s,thresHorz_ft,thresVert_ft,outDir,varargin{:});

%% Assertions
assert(mod(p.Results.encTime_s, p.Results.timeStep_s)==0,'createEncounteres_2:timeMultiples','encTime_s is not a multiple of timeStep_s');
assert(~all(p.Results.isSampleAlt1 && strcmp(p.Results.acmodel1,'bayes')),'createEncounters_2:isSampleAlt1-acmodel1','isSampleAlt1 cannot be true when using bayes model');
assert(~all(p.Results.isSampleAlt2 && strcmp(p.Results.acmodel2,'bayes')),'createEncounters_2:isSampleAlt2-acmodel2','isSampleAlt2 cannot be true when using bayes model');
assert(~all(p.Results.isSampleAlt2 && strcmp(p.Results.acmodel2,'bayes')),'createEncounters_2:isSampleAlt2-acmodel2','isSampleAlt2 cannot be true when using bayes model');
assert(strcmp(p.Results.spheroid.LengthUnit,'foot'),'createEncounters_2:LengthUnit','spheroid LengthUnit must be ''foot''');

%% Initialize encounter structs and encounter counter
waypoints_struct = [];
encCount = 0;

%% Load output from RUN_1
if ~isfile(inFile)
    anchors = table;
    isEncounter = false(0,1);
    files1 = strings(0,1);
    files2 = strings(0,1);
    warning('File not found, inFile = %s...RETURNING',inFile);
    return;
end

load(inFile,'anchors','anchorRange_nm');

% Remove anchors without geospaital trajectories
anchors(anchors.num_geospatial == 0,:) = [];

% Filter by airspace
isAirspace = false(size(anchors,1),1);
for i=1:1:numel(p.Results.classInclude)
    isAirspace = isAirspace | anchors.class == p.Results.classInclude{i};
end
anchors = anchors(isAirspace,:);

if p.Results.anchorPercent ~= 1.0
    pidx = randperm(size(anchors,1), floor(size(anchors,1) * p.Results.anchorPercent));
    anchors = anchors(pidx',:);
end

% Number of anchors
numAnchors = size(anchors,1);

% Tracker if an encounter was generated for anchor
isEncounter = false(numAnchors,1);

% Convert anchorRange_nm from RUN_1 to feet to match spheroid
% An assert above ensures p.Results.spheroid.LengthUnit = 'foot'
anchorRange_ft = anchorRange_nm * unitsratio(p.Results.spheroid.LengthUnit,'nm');

% Set random seed
rng(numAnchors);

%% Return if no data
if numAnchors == 0
    files1 = strings(0,1);
    files2 = strings(0,1);
    fprintf('No anchors in %s, RETURNING\n',inFile);
    return;
end

%% Create clusters
% Clusters currently only used when loading Bayes tracks because
% placeTrack() uses msl2agl(). Loading the DEM for each call of
% placeTrack() is horribly inefficient. To keep DEM size reasonable, we
% break the anchors up into clusters
area_deg = abs(max(anchors.lon_deg) - min(anchors.lon_deg)) * abs(max(anchors.lat_deg) - min(anchors.lat_deg));
area_nm = deg2nm(area_deg);
numClust = ceil(area_nm / 25);
if numel(anchors.lon_deg) <= numClust | numClust == 0
    numClust = ceil(numel(anchors.lon_deg) / 5);
end
[anchors.cluster,~,~] = kmeans([anchors.lon_deg, anchors.lat_deg],numClust);

anchors = sortrows(anchors,{'cluster','num_geospatial','class','lat_deg'});

% Display to screen
fprintf('Starting with %i anchor points across %i clusters\n',numAnchors,numClust);

%% Load FAA DOF
% Load FAA DOF processing in em-core
load(p.Results.dofFile,'Tdof');

% Filter DOF based on obstacle
l = contains(Tdof.obs_type,p.Results.dofObs) & Tdof.alt_ft_agl >= p.Results.dofMinHeight_ft & strcmpi(Tdof.verification_status,'verified');
Tdof = Tdof(l,[4,5,6]); % lat, lon, alt agl;

% Filter DOF based on bounding box of anchor points
bbox = [min(anchors.lon_deg), min(anchors.lat_deg); max(anchors.lon_deg), max(anchors.lat_deg)];
[~, ~, isbox] = filterboundingbox(Tdof.lat_deg,Tdof.lon_deg,bbox);
Tdof = Tdof(isbox,:);

% Create obstacles polygons
[latObstacle,lonObstacle] = scircle1(Tdof.lat_deg,Tdof.lon_deg,repmat(p.Results.dofMaxRange_ft,size(Tdof,1),1),[],p.Results.spheroid,'degrees',20);
altObstacle_ft_agl = Tdof.alt_ft_agl + p.Results.dofMaxVert_ft;

%% Load AC1
switch p.Results.acmodel1
    case 'bayes'
        % Get list of filenames
        listing1 = dir([p.Results.inDir1{1} filesep '*.csv']);
        
        % Concat listings if there is more than one input directory
        if numel(p.Results.inDir2) > 1
            for i=2:1:numel(p.Results.inDir2)
                listing1 = [listing1; dir([p.Results.inDir2{i} filesep '*.csv'])]; %#ok
            end
        end
        assert(~isempty(listing1),'createEncounters_2:listing1','No files found p.Results.inDir2\n');
        % Create full filenames for trajectories which are close to anchor point
        files1 = string(compose('%s/%s',string({listing1.folder}'), string({listing1.name}')));
    case 'geospatial'
        % There could be multiple files for each anchor point
        files1 = strings(sum(anchors.num_geospatial),1);
        s = 1;
        for i=1:1:numAnchors
            e = anchors.num_geospatial(i) + s - 1;
            files1(s:e) = anchors.files{i};
            s = e + 1;;
        end
end

%% Load AC2
switch p.Results.acmodel2
    case 'bayes'
        % Get list of filenames
        listing2 = dir([p.Results.inDir2{1} filesep '*.csv']);
        
        % Concat listings if there is more than one input directory
        if numel(p.Results.inDir2) > 1
            for i=2:1:numel(p.Results.inDir2)
                listing2 = [listing2; dir([p.Results.inDir2{i} filesep '*.csv'])]; %#ok
            end
        end
        assert(~isempty(listing2),'createEncounters_2:listing2','No files found p.Results.inDir2\n');
        % Create full filenames for trajectories which are close to anchor point
        files2 = string(compose('%s/%s',string({listing2.folder}'), string({listing2.name}')));
    case 'geospatial'
        % There could be multiple files for each anchor point
        files2 = strings(sum(anchors.num_geospatial),1);
        s = 1;
        for i=1:1:numAnchors
            e = anchors.num_geospatial(i) + s - 1;
            files2(s:e) = anchors.files{i};
            s = e + 1;;
        end
end

%% Randomly sample a Bayes trajectory for each anchor point
% AC1 must be geospatial
if numel(files1) >= numel(files2)
    % More files in AC1
    files2 = files2(randi(numel(files2),numel(files1),1));
else
    % More files in AC2
    files2 = files2(randperm(numel(files2),numel(files1)));
end

%% Parse out anchor points for each files1/files2 combination
% AC1 must be geospatial
lat0 = zeros(numel(files1),1);
lon0 = zeros(numel(files1),1);
clust0 = zeros(numel(files1),1);
s = 1;
for i=1:1:numAnchors
    e = anchors.num_geospatial(i) + s - 1;
    lat0(s:e) = anchors.lat_deg(i);
    lon0(s:e) = anchors.lon_deg(i);
    clust0(s:e) = anchors.cluster(i);
    s = e + 1;
end

% Update
numAnchors = numel(lat0);

%% Iterate over clusters
for k=1:1:numClust
    % Logical index of cluster
    lk = find(clust0 == k);
    nk = numel(lk);
    
    %% Load Bayes trajectories that are in relative feet
    % Preallocate
    time1_s = cell(size(lk));
    lat1_deg = cell(size(lk));
    lon1_deg = cell(size(lk));
    alt1_ft_msl = cell(size(lk));
    alt1_ft_agl = cell(size(lk));
    isObstacle1 = false(size(lk));
    isShort1 = false(nk,1);
    
    time2_s = cell(size(lk));
    lat2_deg = cell(size(lk));
    lon2_deg = cell(size(lk));
    alt2_ft_msl = cell(size(lk));
    alt2_ft_agl = cell(size(lk));
    isObstacle2 = false(size(lk));
    isShort2 = false(nk,1);
    
    % Aicraft 1
    switch p.Results.acmodel1
        case 'geospatial'
            parfor i=1:1:nk
                % Load trajectory
                [time1_s{i}, lat1_deg{i}, lon1_deg{i}, alt1_ft_msl{i}, alt1_ft_agl{i}] = fastloadtraj(files1(lk(i)));
            end
    end
    
    % Aircraft 2
    switch p.Results.acmodel2
        case 'bayes'
            % Load DEM
            [~,~,Z_m,refvec] = msl2agl(lat0(lk), lon0(lk), p.Results.demBayes,'buff_deg',nm2deg(10));
            if isempty(Z_m)
                [~,~,Z_m,refvec] = msl2agl(lat0(lk), lon0(lk), p.Results.demBayesBackup,'buff_deg',nm2deg(10));
            end
            
            % Iterate over files
            parfor i=1:1:nk
                % Set random seed
                rng(i);
                [time2_s{i},lat2_deg{i},lon2_deg{i},alt2_ft_msl{i},alt2_ft_agl{i},~,~,~,isObstacle2(i)] = placeTrack(files2{lk(i)},lat0(lk(i)), lon0(lk(i)),'z_agl_tol_ft',p.Results.z_agl_tol_ft,...
                    'maxTries',10,...
                    'latObstacle',latObstacle,'lonObstacle',lonObstacle,'altObstacle_ft_agl',altObstacle_ft_agl,...
                    'labelX',p.Results.labelX,'labelY',p.Results.labelY,'labelZ',p.Results.labelZ,...
                    'spheroid',p.Results.spheroid,'dem',p.Results.demBayes,'Z_m',Z_m,'refvec',refvec,'isPlot',false);
            end
            
            isShort2 = cellfun(@numel,time2_s) < p.Results.encTime_s;
        case 'geospatial'
            parfor i=1:1:nk
                % Load trajectory
                [time2_s{i}, lat2_deg{i}, lon2_deg{i}, alt2_ft_msl{i}, alt2_ft_agl{i}] = fastloadtraj(files2(lk(i)));
            end
    end
    %fprintf('Loaded %i pairs for %i / %i clusters\n',nk,k,numClust);
    
    %% Filter if we couldn't avoid obstacles
    isRemove = isObstacle1 | isObstacle2 | isShort1 | isShort2;
    %fprintf('Removing %i pairs due to obstacles or track length\n',sum(isRemove));
    
    time1_s(isRemove,:) = [];
    time2_s(isRemove,:) = [];
    lat1_deg(isRemove,:) = [];
    lat2_deg(isRemove,:) = [];
    lon1_deg(isRemove,:) = [];
    lon2_deg(isRemove,:) = [];
    alt1_ft_msl(isRemove,:) = [];
    alt2_ft_msl(isRemove,:) = [];
    alt1_ft_agl(isRemove,:) = [];
    alt2_ft_agl(isRemove,:) = [];
    
    % Update
    lk(isRemove) = [];
    nk = numel(lk);
    
    % Display status
    fprintf('Valid %i pairs for %i / %i clusters\n',nk,k,numClust);
    
    %% Iterate over unique AC1 / AC2 pairs
    for i=1:1:nk
        % Set random seed
        rng(lk(i));
        
        % Identify which combinaton of waypoints satisfy thresHorz_ft
        [wypts,ac1Time_s,ac1Lat_deg,ac1Lon_deg,ac1Alt_ft_msl,ac2Time_s,ac2Lat_deg,ac2Lon_deg,ac2Alt_ft_msl] = findconflict(p,time1_s{i},lat1_deg{i},lon1_deg{i},alt1_ft_msl{i},time2_s{i},lat2_deg{i},lon2_deg{i},alt2_ft_msl{i},lat0(lk(i)), lon0(lk(i)),anchorRange_ft);
        
        % Check for valid combinations and skip via continue if there are none
        if isempty(wypts)
            %  fprintf('i=%i...No encounter combinations satisfy thresHorz_ft <= %0.1f...skipping via CONTINUE\n',i,thresHorz_ft)
            continue;
        end
        
        % h0_ft origin
        % h0_ft used by geodetic2ned() but we don't use the NED height
        h0_ft = min([ac1Alt_ft_msl;ac2Alt_ft_msl]);
        
        % Sample airspeed and altitude
        ac1Alt_ft_agl = alt1_ft_agl{i}(p.Results.rowstart1:end-p.Results.rowend1); %#ok
        ac2Alt_ft_agl = alt2_ft_agl{i}(p.Results.rowstart2:end-p.Results.rowend2); %#ok
        [v1_kts,v2_kts,altAdjust1_ft,altAdjust2_ft] = samplespeedalt(p,size(wypts,1),ac1Alt_ft_agl,ac2Alt_ft_agl);
        
        % Update waypoints with altitude adjust
        wypts(:,5) = wypts(:,5) + altAdjust1_ft; % AC1 Altitude
        wypts(:,6) = wypts(:,6) + altAdjust2_ft; % AC2 Altitude
        
        % Calculate vertical miss distance for each waypoint
        % Remove waypoints if they violate VMD condition
        vmd = abs(wypts(:,5) - wypts(:,6));
        lvmd = vmd <= p.Results.thresVert_ft;
        wypts = wypts(lvmd,:);
        vmd = vmd(lvmd);
        altAdjust1_ft = altAdjust1_ft(lvmd);
        altAdjust2_ft = altAdjust2_ft(lvmd);
        
        % Check for valid combinations and skip via continue if there are none
        if isempty(wypts)
            continue;
        end
        
        % Convert each aircraft to north / east coordinates in feet
        % The CSIM waypoint uses feet but code mostly calculates using lat / lon
        % We're currently ignoring the altitude adjustment because this
        % function is being used for low altitudes...calling this once is
        % faster than putting it in the j for loop below
        [ac1North_ft,ac1East_ft,~] = geodetic2ned(ac1Lat_deg,ac1Lon_deg,ac1Alt_ft_msl,lat0(lk(i)),lon0(lk(i)),h0_ft,p.Results.spheroid);
        [ac2North_ft,ac2East_ft,~] = geodetic2ned(ac2Lat_deg,ac2Lon_deg,ac2Alt_ft_msl,lat0(lk(i)),lon0(lk(i)),h0_ft,p.Results.spheroid);
        
        % Check to make sure coordinates are not NaN
        if any(isnan([ac1North_ft; ac1North_ft; ac2North_ft; ac2East_ft]))
            warning('createEncounters_2:geodetic2ned','geodetic2ned() outputs NaN coordinates\nTrajectory files:\nAC1 = %s\nAC2 = %s\nSkipping via CONTINUE\n',files1(lk(i)),files2(lk(i)));
            continue;
        end
        
        % Iterate over unique waypoint pairs
        encAdd = false;
        jadd = false(size(wypts,1),1);
        for j=1:1:size(wypts,1)
            % Crop tracks
            [ac1NEU_foward,ac2NEU_foward,ac1NEU_backward,ac2NEU_backward] = croptracks(p,wypts(j,:),encTime_s,ac1Time_s,ac2Time_s,ac1Lat_deg,ac2Lat_deg,ac1Lon_deg,ac2Lon_deg,ac1North_ft,ac2North_ft,ac1East_ft,ac2East_ft,ac1Alt_ft_msl,ac2Alt_ft_msl,altAdjust1_ft(j),altAdjust2_ft(j),v1_kts(j),v2_kts(j));
            
            % Update waypoints_csim struct
            [waypoints_struct,encCount,encAdd] = updatewaypointstruct(p,waypoints_struct,encCount,ac1NEU_foward,ac2NEU_foward,ac1NEU_backward,ac2NEU_backward);
            jadd(j) = encAdd;
        end % End j loop
        
        % Record if an encounter was made for anchor point
        if encAdd
            isEncounter(lk(i)) = true;
            wypts = wypts(jadd,:);
        end

        if p.Results.isPlot && size(wypts,1) > 0
            if ~encAdd
                %close(lk(i));
            else
                % Create figure and basemap
                figure(lk(i));
                gx = geoaxes('Basemap','topographic',...
                    'FontSize',16,...
                    'FontWeight','bold',...
                    'Position',[0 0 1 1]);
                bmap = geobasemap(gx);
                geolimits(lat0(lk(i)) + nm2deg([-anchorRange_nm anchorRange_nm]*1.05), lon0(lk(i)) + nm2deg([-anchorRange_nm anchorRange_nm]*1.05));
                hold on;
                
                % Anchor Radius
                [lat0c,lon0c] = scircle1(lat0(lk(i)),lon0(lk(i)),anchorRange_ft,[],p.Results.spheroid,'degrees',100);
                geoplot(lat0c,lon0c,'k--');
                
                % Paths
                geoplot(ac1Lat_deg,ac1Lon_deg,'b-');
                geoplot(ac2Lat_deg,ac2Lon_deg,'r-');
                
                % Anchor Point
                geoscatter(lat0(lk(i)),lon0(lk(i)),50,'k','d','filled');
                
                % Threshold Points
                for w=1:1:size(wypts,1)
                    if w==1; hv = 'on'; else hv = 'off'; end; % Only want to plot conflict points once on legend
                    geoscatter(wypts(w,1),wypts(w,2),60,'b','h','filled','HandleVisibility',hv)
                    geoscatter(wypts(w,3),wypts(w,4),60,'r','p','filled','HandleVisibility',hv)
                    
                    % Alternate horizontal alignment for odd / even values of w
                    % o1, o2 = offsets for text
                    if mod(w,2) == 0
                        ha1 = 'left';
                        ha2 = 'right';
                        o1 = nm2deg(0.04);
                        o2 = -nm2deg(0.04);
                    else
                        ha1 = 'right';
                        ha2 = 'left';
                        o1 = -nm2deg(0.04);
                        o2 = nm2deg(0.04);
                    end
                    
                    text(wypts(w,1),wypts(w,2)+o1,num2str(w),'Color','b','FontSize',8,'FontWeight','bold','HorizontalAlignment',ha1);
                    text(wypts(w,3),wypts(w,4)+o2,num2str(w),'Color','r','FontSize',8,'FontWeight','bold','HorizontalAlignment',ha2);
                end
                
                % Adjust map size, legend, hold off
                set(gcf,'Units','inches','Position',[1 1 11.94 5.28])
                legend(sprintf('Radius = %0.2fnm',anchorRange_nm),sprintf('AC1 = %s',p.Results.acmodel1), sprintf('AC2 = %s',p.Results.acmodel2),sprintf('(%0.3f,%0.3f)',lat0(lk(i)),lon0(lk(i))),'AC1 @ conflict','AC2 @ conflict');
                hold off;
                
                % Print
                print(['output' filesep 'plot_encounter_examples' filesep 'example_2' '_lk' num2str(lk(i)) '_k' num2str(k) '_i' num2str(i) '_lat' num2str(lat0(lk(i))) '_lon' num2str(lon0(lk(i))) '_thresHorz' num2str(thresHorz_ft) '_thresVert' num2str(thresVert_ft) '_ac1' p.Results.acmodel1 '_ac2'  p.Results.acmodel2 '.png'],'-dpng','-r300');
            end
        end
    end % End i loop
    fprintf('%i encounters after %i / %i clusters\n',encCount,k,numClust);
end % End k loop

%% Write waypoints to file
% Make output directory
mkdir(outDir);

% Write individual encounters to file
parfor files_idx = 1:1:size(waypoints_struct,2)
    output_enc_filename = [outDir filesep 'enc_' num2str(files_idx) '.dat'];
    save_waypoints(output_enc_filename,waypoints_struct(:,files_idx));
end

% Display status to screen
fprintf('%i encounters written to file\n',size(waypoints_struct,2));

%% Zip waypoints for easy sharing
if encCount > 0
    if p.Results.isZip; zip([outDir '.zip'],[outDir filesep '*.dat']); end
else
    fprintf('Nothing to zip, no encounter generated for %s\n',inFile);
end
