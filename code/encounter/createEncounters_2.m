function [metadata] = createEncounters_2(inAnchors,outDir,encTime_s,thresHorz_ft,thresVert_ft,varargin)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
% CREATEENCOUNTERS_2 OUTPUTS ENCOUNTERS IN WAYPOINT FORMAT
%
% SEE ALSO loadTrack sampleSpeedAlt adjustTrack findConflict findCropIdx trackTimetable2NEU plotEncounter

%% Input parser
p = inputParser;

% Required
addRequired(p,'inAnchors'); % Output of RUN_1
addRequired(p,'outDir'); % Output directory
addRequired(p,'encTime_s',@isnumeric); % Duration of encounter
addRequired(p,'thresHorz_ft',@isnumeric);
addRequired(p,'thresVert_ft',@isnumeric);

% Optional - Minimum encounter initial conditions
addParameter(p,'conflictTime_s',ceil(encTime_s/2),@isnumeric); % Time at which conflict critieria are satisfied
addParameter(p,'initHorz_ft',[6076 97217],@(x) isnumeric(x) && numel(x) == 2); %97217 ~ 16 nm ~ 421 feet (250 knots) * 240 seconds (uncor MVL)
addParameter(p,'initVert_ft',[0 1500],@(x) isnumeric(x) && numel(x) == 2);
addParameter(p,'spheroid_ft',wgs84Ellipsoid('ft'));

% Optional - Anchor
addParameter(p,'classInclude',{'B','C','D','O'},@iscell);

% Optional - Encounter Variations
addParameter(p,'maxTries',6,@(x) isnumeric(x) && numel(x) == 1);
addParameter(p,'maxEncPerPair',1,@isnumeric); % Maximum waypoint pairs for each unique pair
addParameter(p,'timeStep',seconds(1),@isduration); %  Output timestep, CSIM requires at least >= 1
addParameter(p,'anchorPercent',1.0,@isnumeric); %  Maximum waypoint combinations to consider for each pair

% Optional - AC permutations
addParameter(p,'mdlType1','geospatial',@(x) ischar(x) && any(strcmpi(x,{'bayes','geospatial'})));
addParameter(p,'mdlType2','bayes',@(x) ischar(x) && any(strcmpi(x,{'bayes','geospatial'})));
addParameter(p,'inDir1',{},@iscell);
addParameter(p,'inDir2',{},@iscell);

% Optional - Airspeed Sampling (Geospatial only)
addParameter(p,'isSampleV1',true,@islogical);
addParameter(p,'isSampleV2',false,@islogical);
addParameter(p,'rangeV1_ft_s',[9 168],@(x) isnumeric(x) && numel(x) == 2); % ~ [5,100] knots
addParameter(p,'rangeV2_ft_s',[9 422],@(x) isnumeric(x) && numel(x) == 2); % ~ [5,100] knots

% Optional - Altitude Sampling (Geospatial only)
addParameter(p,'isSampleAlt1',true,@islogical);
addParameter(p,'isSampleAlt2',false,@islogical);
addParameter(p,'rangeAlt1_ft_agl',[50 1200],@(x) isnumeric(x) && numel(x) == 2);
addParameter(p,'rangeAlt2_ft_agl',[50 2000],@(x) isnumeric(x) && numel(x) == 2);
addParameter(p,'isGeoPointOverride',true,@islogical);

% Optional - Bayes Tracks
addParameter(p,'bayesFile1',[getenv('AEM_DIR_BAYES') filesep 'model' filesep 'uncor_1200only_fwse_v1p2.txt']);
addParameter(p,'bayesFile2',[getenv('AEM_DIR_BAYES') filesep 'model' filesep 'uncor_1200only_fwse_v1p2.txt']);
addParameter(p,'z_agl_tol_ft',200,@(x) isnumeric(x) && numel(x) == 1); % Used by placeTrack()
addParameter(p,'dem','dted1',@(x) isstr(x) && any(strcmpi(x,{'dted1','dted2','globe','gtopo30','srtm1','srtm3','srtm30'})));
addParameter(p,'demDir',char.empty(0,0),@ischar);
addParameter(p,'demBackup','globe',@(x) isstr(x) && any(strcmpi(x,{'dted1','dted2','globe','gtopo30','srtm1','srtm3','srtm30'})));
addParameter(p,'demBackupDir',char.empty(0,0),@ischar);
addParameter(p,'labelTime','time_s');
addParameter(p,'labelX','x_ft');
addParameter(p,'labelY','y_ft');
addParameter(p,'labelZ','z_ft');

% Optional - FAA DOF
addParameter(p,'dofFile',[getenv('AEM_DIR_CORE') filesep 'output' filesep 'dof.mat'],@ischar); % Location of output of RUN_readfaadof from em-core
addParameter(p,'dofTypes',{'ag equip','antenna','lgthouse','met','monument','silo','spire','stack','t-l twr','tank','tower','utility pole','windmill'},@iscell); % Obstacles to keep
addParameter(p,'dofMinHeight_ft',50,@isnumeric); % Minimum height of obstacles
addParameter(p,'dofMaxRange_ft',500,@isnumeric);
addParameter(p,'dofMaxVert_ft',50,@isnumeric);

% Optional - Random seed
addParameter(p,'initialSeed',0,@isnumeric);

% Optional - Output
addParameter(p,'outHash',char.empty(0,0),@ischar);
addParameter(p,'isPlot',false,@islogical);
addParameter(p,'isZip',false,@islogical);

% Row offsets
% These are needed because the geospatial trajectory generator will add
% vertical take off and landing segements, which you may not want
addParameter(p,'isStartEndNonZero1',true,@islogical); %
addParameter(p,'isStartEndNonZero2',true,@islogical); %

% Parse
parse(p,inAnchors,outDir,encTime_s,thresHorz_ft,thresVert_ft,varargin{:});

initHorz_ft = p.Results.initHorz_ft;
initVert_ft = p.Results.initVert_ft;
conflictTime_s = p.Results.conflictTime_s;

dofMaxRange_ft = p.Results.dofMaxRange_ft;
dofMaxVert_ft = p.Results.dofMaxVert_ft;

maxTries = p.Results.maxTries;
maxEncPerPair = p.Results.maxEncPerPair;
timeStep = p.Results.timeStep;
spheroid_ft = p.Results.spheroid_ft;
z_agl_tol_ft = p.Results.z_agl_tol_ft;
mdlType1 = p.Results.mdlType1;
mdlType2 = p.Results.mdlType2;
isPlot = p.Results.isPlot;

isStartEndNonZero1 = p.Results.isStartEndNonZero1;
isStartEndNonZero2 = p.Results.isStartEndNonZero2;

isSampleAlt1 = p.Results.isSampleAlt1;
isSampleV1 = p.Results.isSampleV1;
isSampleAlt2 = p.Results.isSampleAlt2;
isSampleV2 = p.Results.isSampleV2;
isGeoPointOverride = p.Results.isGeoPointOverride;

rangeAlt1_ft_agl = p.Results.rangeAlt1_ft_agl;
rangeV1_ft_s = p.Results.rangeV1_ft_s;
rangeAlt2_ft_agl = p.Results.rangeAlt2_ft_agl;
rangeV2_ft_s = p.Results.rangeV2_ft_s;

%% Assertions and More Input Handling
if strcmpi(mdlType1,'bayes') & (isSampleAlt1 | isSampleV1)
    isSampleAlt1 = false;
    isSampleV1 = false;
    warning('Cannot adjust altitude and velocity for Bayes tracks, setting isSampleAlt1 and isSampleV1 to false');
end

if strcmpi(mdlType2,'bayes') & (isSampleAlt2 | isSampleV2)
    isSampleAlt2 = false;
    isSampleV2 = false;
    warning('Cannot adjust altitude and velocity for Bayes tracks, setting isSampleAlt2 and isSampleV2 to false');
end

assert(all(encTime_s >= conflictTime_s),'Failed to satisfy encTime_s >= conflictTime_s');
assert(~all(isSampleAlt1 && strcmp(mdlType1,'bayes')),'createEncounters_2:isSampleAlt1-acmodel1','isSampleAlt1 cannot be true when using bayes model');
assert(~all(isSampleAlt2 && strcmp(mdlType2,'bayes')),'createEncounters_2:isSampleAlt2-acmodel2','isSampleAlt2 cannot be true when using bayes model');
assert(strcmp(spheroid_ft.LengthUnit,'foot'),'createEncounters_2:LengthUnit','spheroid LengthUnit must be ''foot''');

% Look ahead and back time relative to conflict time
durBeforeConflict_s = conflictTime_s;
durAfterConflict_s = encTime_s - conflictTime_s;

%% Inputs hardcode
bayesDuration_s = 240;
maxSeed = 2^32;
encPerFile = 100000;
metaVarNames = {'encId','anchorRow','anchorFile','lat0_deg','lon0_deg','el0_ft_msl','randSeed',...
    'hmd_conflict_ft','vmd_conflict_ft','hmd_cpa_ft','vmd_cpa_ft','hmd_initial_ft','vmd_initial_ft','hmd_final_ft','vmd_final_ft',...
    'speed1_conflict_ft_s','speed2_conflict_ft_s','speed1_cpa_ft_s','speed2_cpa_ft_s','speed1_initial_ft_s','speed2_initial_ft_s','speed1_final_ft_s','speed2_final_ft_s'...
    'time_enc_s','time_conflict_s','time_cpa_s',...
    'G','A'};

%% Initialize encounter structs and encounter counter
encCount = 0;
randSeed = p.Results.initialSeed;

if isempty(p.Results.outHash)
    outName = 'enc_';
else
    outName = ['enc_' p.Results.outHash '_'];
end

rng(randSeed,'twister');

%% Load output from RUN_1

if ischar(inAnchors) && isfile(inAnchors)
    
    % Load anchors
    load(inAnchors,'anchors','anchorRange_nm');
    
    % Remove anchors without geospaital trajectories
    anchors(anchors.num_geospatial == 0,:) = [];
    
    % Filter by airspace class
    isAirspace = false(size(anchors,1),1);
    for i=1:1:numel(p.Results.classInclude)
        isAirspace = isAirspace | anchors.class == p.Results.classInclude{i};
    end
    anchors = anchors(isAirspace,:);
    
    % If both aircraft are geospatial, filter to anchors with multiple files
    if strcmpi(mdlType1,'geospatial') && strcmpi(mdlType2,'geospatial')
        anchors =  anchors(anchors.num_geospatial > 1,:);
    end
    
    % Filter by anchor percentage
    if p.Results.anchorPercent ~= 1.0
        pidx = randperm(size(anchors,1), floor(size(anchors,1) * p.Results.anchorPercent));
        anchors = anchors(pidx',:);
    end
    
else
    if istable(inAnchors)
        anchors = inAnchors;
        anchorRange_nm = inAnchors.anchorRange_nm(1);
    else
        metadata = table.empty(0,numel(metaVarNames));
        warning('Input anchors not found...RETURNING');
        return;
    end
end

% Number of anchors
numAnchors = size(anchors,1);

% Convert anchorRange_nm from RUN_1 to feet to match spheroid
% An assert above ensures p.Results.spheroid.LengthUnit = 'foot'
anchorRange_ft = anchorRange_nm * unitsratio(spheroid_ft.LengthUnit,'nm');

%% Update anchors with geographic domain and airspace class
% This information is used to set the UncorEncounterModel start property
if strcmpi(mdlType1,'bayes') | strcmpi(mdlType2,'bayes')
    % Geographic domain
    if ~any(strcmpi(anchors.Properties.VariableNames,'G'))
        warning('anchors table does not have variable, G; calculating now');
        %anchors.G = ones(size(anchors,1),1);
        anchors.G = IdentifyGeographicVariable(anchors.lat_deg,anchors.lon_deg);
    end
    
    % Airspace class
    anchors.A = repmat(4,size(anchors,1),1);
    anchors.A(anchors.class == 'B') = 1;
    anchors.A(anchors.class == 'C') = 2;
    anchors.A(anchors.class == 'D') = 3;
    anchors.A(anchors.class == 'O') = 4;
end

%% Load FAA DOF
BoundingBox_wgs84 = [min(anchors.lon_deg), min(anchors.lat_deg); max(anchors.lon_deg), max(anchors.lat_deg)];

[~, Tdof] = gridDOF('inFile',p.Results.dofFile,...
    'BoundingBox_wgs84',BoundingBox_wgs84,...
    'minHeight_ft',p.Results.dofMinHeight_ft,...
    'isVerified',true,...
    'obsTypes',p.Results.dofTypes);

%% Load Bayes models
if strcmpi(mdlType1,'bayes'); mdl1 = UncorEncounterModel('parameters_filename',p.Results.bayesFile1); end
if strcmpi(mdlType2,'bayes'); mdl2 = UncorEncounterModel('parameters_filename',p.Results.bayesFile2); end

%% Preallocate metadata
nEnc = sum(anchors.num_geospatial) * maxEncPerPair;
pad = nan(nEnc,1); %pad = PreAllocated Double
metadata = table((1:1:nEnc)',pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,repmat(encTime_s,nEnc,1),repmat(conflictTime_s,nEnc,1),pad,pad,pad,...
    'VariableNames',metaVarNames);

%% Iterate over anchors
cEnc = 0; % Encounter counter
cFiles = 0; % Files written counter
waypointsBuffer = [];

for ii=1:1:numAnchors
    % Parse iith anchor
    iiAnchor = anchors(ii,:);
    lat0_deg = iiAnchor.lat_deg;
    lon0_deg = iiAnchor.lon_deg;
    iiFiles = iiAnchor.files{1};
    
    % Load DEM near anchor
    % 1 knot = 1.68781 feet per second
    rad_ft = encTime_s * 600 * 1.68781;
    [latc, lonc] = scircle1(lat0_deg,lon0_deg,rad_ft,[],spheroid_ft);
    dem = p.Results.dem;
    demDir = p.Results.demDir;
    [el0_ft_msl,~,Z_m,refvec,R] = msl2agl([lat0_deg; min(latc); max(latc)], [lon0_deg; min(lonc); max(lonc)],dem,'demDir',demDir,'isCheckOcean',true,'isFillAverage',true);
    
    % If primary DEM didn't work, try backup
    % This is often due to missing DEM tiles
    if isempty(el0_ft_msl)
        dem = p.Results.demBackup;
        demDir = p.Results.demBackupDir;
        [el0_ft_msl,~,Z_m,refvec,R] = msl2agl([lat0_deg; min(latc); max(latc)], [lon0_deg; min(lonc); max(lonc)],dem,'demDir',demDir,'isCheckOcean',true,'isFillAverage',true);
    end
    
    % Get MSL elevation of (lat0_deg,lon0_deg)
    el0_ft_msl = el0_ft_msl(1);
    
    % Filter obstacles within range of anchor
    bbox = [min(lonc), min(latc); max(lonc), max(latc)];
    [~, ~, isbox] = filterboundingbox(Tdof.lat_deg,Tdof.lon_deg,bbox);
    iiDof = Tdof(isbox,:);
    
    %     if isPlot
    %          fig = plotAnchor(iiAnchor,iiDof);
    %     end
    
    % Start distribution for sampling uncorrelated encounter model
    if strcmpi(mdlType1,'bayes') | strcmpi(mdlType2,'bayes')
        bayesStart = {iiAnchor.G,iiAnchor.A,[],[],[],[],[]};
    end
    
    % Altitude start distribution is random and set further in the code
    if strcmpi(mdlType1,'bayes');
        idxL1 = find(strcmp(mdl1.labels_initial,'"L"'));
        binMaxL1 = find(mdl1.boundaries{idxL1} < rangeAlt1_ft_agl(2),1,'last');
    end
    
    if strcmpi(mdlType2,'bayes');
        idxL2 = find(strcmp(mdl2.labels_initial,'"L"'));
        binMaxL2 = find(mdl2.boundaries{idxL2} < rangeAlt2_ft_agl(2),1,'last');
    end
    
    % Create pairs
    if strcmpi(mdlType1,'geospatial') && strcmpi(mdlType2,'geospatial')
        pairs =  nchoosek(1:iiAnchor.num_geospatial,2);
    else
        % both bayes
        if strcmpi(mdlType1,'bayes') && strcmpi(mdlType2,'bayes')
            pairs = [nan, nan];
        else
            % mdlType1 = 'geospatial' and mdlType2 = 'bayes'
            if strcmpi(mdlType1,'geospatial')
                pairs = [(1:iiAnchor.num_geospatial)', nan(iiAnchor.num_geospatial,1)];
            else
                % mdlType1 = 'bayes' and mdlType2 = 'geospatial'
                pairs = [nan(iiAnchor.num_geospatial,1), (1:iiAnchor.num_geospatial)'];
            end
        end
    end
    
    % Iterate over pairs
    for jj=1:1:size(pairs,1)
        
        % Load Aircraft 1
        switch mdlType1
            case 'geospatial'
                track1 = loadTrack(mdlType1, iiFiles(pairs(jj,1)), isStartEndNonZero1,rangeAlt1_ft_agl, rangeV1_ft_s, spheroid_ft, lat0_deg,lon0_deg,el0_ft_msl);
            case 'bayes'
                % Update Bayesian start distribution
                bayesStart{idxL1} = randi(binMaxL1,1);
                mdl1.start = bayesStart;
                track1 = loadTrack(mdlType1, mdl1, isStartEndNonZero1, rangeAlt1_ft_agl, rangeV1_ft_s, spheroid_ft, lat0_deg,lon0_deg,el0_ft_msl,Z_m,refvec,R,iiDof,dofMaxRange_ft,dofMaxVert_ft,bayesDuration_s,z_agl_tol_ft);
        end
        
        % Initialize counters
        c = 1;
        jjNumEncs = 0;
        conflicts = double.empty(0,10);
        
        % Iterate until we get the desired number of encounters or we reach  maxTries
        while c <= maxTries
            
            % Reload Aircraft 1 if Bayes
            switch mdlType1
                case 'bayes'
                    % Update Bayesian start distribution
                    bayesStart{idxL1} = randi(binMaxL1,1);
                    mdl1.start = bayesStart;
                    track1 = loadTrack(mdlType1, mdl1, isStartEndNonZero1, rangeAlt1_ft_agl, rangeV1_ft_s, spheroid_ft, lat0_deg,lon0_deg,el0_ft_msl,Z_m,refvec,R,iiDof,dofMaxRange_ft,dofMaxVert_ft,bayesDuration_s,z_agl_tol_ft);
            end
            
            % Load Aircraft 2
            switch mdlType2
                case 'geospatial'
                    track2 = loadTrack(mdlType2, iiFiles(pairs(jj,2)), isStartEndNonZero2, rangeAlt2_ft_agl, rangeV2_ft_s, spheroid_ft, lat0_deg,lon0_deg,el0_ft_msl);
                case 'bayes'
                    % Update Bayesian start distribution
                    bayesStart{idxL2} = randi(binMaxL2,1);
                    mdl2.start = bayesStart;
                    track2 = loadTrack(mdlType2, mdl2, isStartEndNonZero2, rangeAlt2_ft_agl, rangeV2_ft_s, spheroid_ft, lat0_deg,lon0_deg,el0_ft_msl,Z_m,refvec,R,iiDof,dofMaxRange_ft,dofMaxVert_ft,bayesDuration_s,z_agl_tol_ft);
            end
            
            % Determine if tracks are short (this is bad)
            isShort1 = size(track1,1) < encTime_s  & ~isSampleV1;
            isShort2 = size(track2,1) < encTime_s  & ~isSampleV2;
            
            % Do something if not empty, not missing, and not short
            if ~isempty(track1) && ~isempty(track2) && ~any(ismissing(track1),'all') && ~any(ismissing(track2),'all') && ~isShort1 && ~isShort2
                
                % Sample airspeed and altitude adjustments
                % If desired, override geospatial tracks assumed to be
                % point inspections. Since sampleSpeedAlt() checks model
                % type, Bayes tracks will still return no adjustements
                pointPatterns = {'uswtdb','dof'};
                if isGeoPointOverride && strcmpi(mdlType1,'geospatial') && contains(iiFiles(pairs(jj,1)),pointPatterns,'IgnoreCase',true)
                    [altAdjust1_ft, vAdjust1_ft_s] = sampleSpeedAlt(mdlType1,track1, rangeAlt1_ft_agl, [2 16], false, true, 1);
                else
                    [altAdjust1_ft, vAdjust1_ft_s] = sampleSpeedAlt(mdlType1,track1, rangeAlt1_ft_agl, rangeV1_ft_s, isSampleAlt1, isSampleV1, 1);
                end
                
                if isGeoPointOverride && strcmpi(mdlType2,'geospatial') && contains(iiFiles(pairs(jj,2)),pointPatterns,'IgnoreCase',true)
                    [altAdjust2_ft, vAdjust2_ft_s] = sampleSpeedAlt(mdlType2,track2, rangeAlt2_ft_agl, [2 16], false, true, 1);
                else
                    [altAdjust2_ft, vAdjust2_ft_s] = sampleSpeedAlt(mdlType2,track2, rangeAlt2_ft_agl, rangeV2_ft_s, isSampleAlt2, isSampleV2, 1);
                end
                
                % Adjust altitude and speed
                ctrack1 = adjustTrack(track1, altAdjust1_ft, vAdjust1_ft_s);
                ctrack2 = adjustTrack(track2, altAdjust2_ft, vAdjust2_ft_s);
                
                % Identify which combinaton of waypoints satisfy
                % HMD, VMD, and time criteria
                conflicts = findConflict(mdlType1,mdlType2,ctrack1,ctrack2,lat0_deg,lon0_deg,anchorRange_ft,thresHorz_ft,thresVert_ft,durBeforeConflict_s,durAfterConflict_s,spheroid_ft);
                
                % Do something if there is any conflicts
                if ~isempty(conflicts)
                    
                    % Find combination of start and end indicies
                    [conflicts,sidx1,sidx2,eidx1,eidx2] = findCropIdx(conflicts,mdlType1,mdlType2,ctrack1,ctrack2,initHorz_ft,initVert_ft,durBeforeConflict_s,durAfterConflict_s);
                    
                    % Do something if conflicts are not empty
                    if ~isempty(conflicts)
                        % Select one
                        idx = randi(size(conflicts,1),1,1);
                        conflicts = conflicts(idx,:);
                        sidx1 = sidx1(idx);
                        sidx2 = sidx2(idx);
                        eidx1 = eidx1(idx);
                        eidx2 = eidx2(idx);
                        
                        % Convert timetables to arrays in a specific column order
                        [ac1NEU, ac1speed_ft_s] = trackTimetable2NEU(ctrack1,sidx1,eidx1,timeStep);
                        [ac2NEU, ac2speed_ft_s] = trackTimetable2NEU(ctrack2,sidx2,eidx2,timeStep);
                        
                        % Calculate horizontal and vertical seperation
                        r_ft = hypot(ac1NEU(:,2)-ac2NEU(:,2),ac1NEU(:,3)-ac2NEU(:,3));
                        z_ft = abs(ac1NEU(:,4)-ac2NEU(:,4));
                        d_ft = sqrt((ac1NEU(:,2)-ac2NEU(:,2)).^2 + (ac1NEU(:,3)-ac2NEU(:,3)).^2 + (ac1NEU(:,4)-ac2NEU(:,4)).^2);
                        
                        % Confirm that initial conditions are satisifed
                        initHMDGood = r_ft(1) >= initHorz_ft(1) & r_ft(1) <= initHorz_ft(2);
                        initVMDGood = z_ft(1) >= initVert_ft(1) & z_ft(1) <= initVert_ft(2);
                        if initHMDGood & initVMDGood
                            
                            % Identify when CPA occurs
                            [~,idxCPA] = min(d_ft);
                            
                            % Create waypoints struct
                            cwaypoints = neu2Waypoints(ac1NEU,ac2NEU);
                            
                            % Append waypoints buffer
                            waypointsBuffer = [waypointsBuffer cwaypoints];
                            
                            % Update counters
                            cEnc = cEnc + 1;
                            jjNumEncs = jjNumEncs + 1;
                            
                            % Update metadata
                            metadata.anchorRow(cEnc) = ii;
                            metadata.anchorFile(cEnc) = jj;
                            metadata.lat0_deg(cEnc) = lat0_deg;
                            metadata.lon0_deg(cEnc) = lon0_deg;
                            metadata.el0_ft_msl(cEnc) = el0_ft_msl;
                            metadata.randSeed(cEnc) = randSeed;
                            
                            % Update metadata: HMD / VMD
                            metadata.hmd_conflict_ft(cEnc) = conflicts(:,9);
                            metadata.vmd_conflict_ft(cEnc) = conflicts(:,10);
                            metadata.hmd_cpa_ft(cEnc) = r_ft(idxCPA);
                            metadata.vmd_cpa_ft(cEnc) = z_ft(idxCPA);
                            metadata.hmd_initial_ft(cEnc) = r_ft(1);
                            metadata.vmd_initial_ft(cEnc) = z_ft(1);
                            metadata.hmd_final_ft(cEnc) = r_ft(end);
                            metadata.vmd_final_ft(cEnc) = z_ft(end);
                            
                            % Update metadata: speed
                            metadata.speed1_conflict_ft_s(cEnc) = ac1speed_ft_s(conflictTime_s);
                            metadata.speed2_conflict_ft_s(cEnc) = ac2speed_ft_s(conflictTime_s);
                            metadata.speed1_cpa_ft_s(cEnc) = ac1speed_ft_s(idxCPA);
                            metadata.speed2_cpa_ft_s(cEnc) = ac2speed_ft_s(idxCPA);
                            metadata.speed1_initial_ft_s(cEnc) = ac1speed_ft_s(1);
                            metadata.speed2_initial_ft_s(cEnc) = ac2speed_ft_s(1);
                            metadata.speed1_final_ft_s(cEnc) = ac1speed_ft_s(end);
                            metadata.speed2_final_ft_s(cEnc) = ac2speed_ft_s(end);
                            
                            % Update metadata: cpa time
                            metadata.time_cpa_s(cEnc) = idxCPA;
                            
                            % Update metadata: G, A
                            metadata.G(cEnc) = iiAnchor.G;
                            metadata.A(cEnc) = iiAnchor.A;
                            
                            % Write to file
                            if cEnc > 0 && mod(cEnc,encPerFile) == 0
                                % Create filename and write to file
                                outFile = [outDir filesep outName sprintf('%05.f',cFiles) '.dat'];
                                save_waypoints(outFile,waypointsBuffer);
                                
                                % Display status to screen
                                fprintf('%i encounters written to file: %s\n',cEnc,outFile);
                                
                                % Update counter and reset buffer
                                cFiles = cFiles + 1;
                                waypointsBuffer = [];
                            end
                            
                            % Plot
                            if isPlot
                                fig = plotEncounter(conflicts,mdlType1,mdlType2,ctrack1,ctrack2,sidx1,eidx1,sidx2,eidx2,lat0_deg,lon0_deg,randSeed);
                            end
                        end
                    end
                end
            end
            
            % Advance counter
            % If sufficient encounters generated, set counter to Inf to break loop
            if jjNumEncs >= maxEncPerPair
                nTries = c;
                c = inf;
            else
                nTries = c;
                c = c + 1;
            end
            
            % Update random seed
            randSeed = randSeed + 1;
            if randSeed < maxSeed
                rng(randSeed,'twister');
            else
                warning('Maximum random seed exceeded, seed not updated');
            end
            
        end
        
        % Display status
        fprintf('%i encounters generated in %i attempts when lat0_deg=%0.4f, lon0_deg=%0.4f, ii=%i, jj=%i\n',jjNumEncs,nTries,lat0_deg,lon0_deg,ii,jj);
    end
    
    % Display status
    %fprintf('%i / %i\n',ii,numAnchors);
end % End ii


% Write last file of remaining encounters in waypointsBuffer
if cEnc > 0 && ~isempty(waypointsBuffer)
    % Create filename and write to file
    outFile = [outDir filesep outName sprintf('%05.f',cFiles) '.dat'];
    save_waypoints(outFile,waypointsBuffer);
    
    % Display status to screen
    fprintf('%i encounters written to file: %s\n',cEnc,outFile);
end

% Remove unused rows in metadata table
if cEnc > 0
    metadata = metadata(1:cEnc,:);
else
    metadata = table.empty(0,13);
end

%% Zip waypoints for easy sharing
if cEnc > 0
    if p.Results.isZip; zip([outDir '.zip'],[outDir filesep '*.dat']); end
else
    fprintf('Nothing to zip, no encounter generated for %s\n',inAnchors);
end
