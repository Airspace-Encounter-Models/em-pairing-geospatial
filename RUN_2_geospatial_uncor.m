% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
% Tradespace
tradespace = CreateTradespacePairingGeo('mdlType1','geospatial','mdlType2','bayes',...
    'thresHorz_ft',4500,'thresVert_ft',450,...
    'dofMaxRange_ft',[500],'dofMaxVert_ft',[500],...
    'conflictTime_s',[60:30:210, 60:30:180]','encTime_s',[(60:30:210)+30, (60:30:180)+60]');
tradespace = sortrows(tradespace,{'dofMaxRange_ft','encTime_s'});
        
% Parent output directory
outDirRoot = ['output' filesep '2_encounters' filesep 'geospatial_uncor' filesep datestr(date,'yyyy-mm-dd')];

% RUN_1 parameters
anchorRange_nm = 1.00;
usecase = 'all';

seedGen = 1;

% Other parameter not defined in tradespace table
dem = 'globe';

%% Suppress warning
warning('off','map:io:UnableToDetermineCoordinateSystemType')

%% Create tasks
[tasks, anchors,initialSeeds] = CreateTasksPairingGeo(tradespace,anchorRange_nm,usecase,seedGen);

%% Save tradespace
[~, ~, ~] = mkdir(outDirRoot);
save([outDirRoot filesep 'tradespace.mat'],'tradespace','tasks');

%% Create status table to record performance
status = table( tradespace.configId,zeros(size(tradespace,1),1),zeros(size(tradespace,1),1),'VariableNames',{'tradespaceId','workTime_s','nEncounters'});

%% Iterate over tasks
for ii=1:1:size(tasks,1)
    % Start timer
    tic
    
    % Display status
    fprintf('Starting task global id %i\n',ii);
    
    % Filter anchors
    iiAnchors = anchors(tasks.sidx(ii):tasks.eidx(ii),:);
    
    iiGlobal = tasks.globalId(ii);
    iiConfig = tasks.configId(ii);
    
    % Tradespace idx
    idxTS = find( tradespace.configId == iiConfig);
    
    % Create output directory
    outHash = sprintf('%04.f_%07.f',iiConfig,iiGlobal);
    outDir = [outDirRoot filesep sprintf('%04.f',iiConfig)];
    [~, ~, ~] = mkdir(outDir);
     
    % Create encounters
    metadata = createEncounters_2(iiAnchors,outDir,...
        tradespace.encTime_s(idxTS,:),tradespace.thresHorz_ft(idxTS,:),tradespace.thresVert_ft(idxTS,:),...
        'bayesFile1',tradespace.bayesFile1{idxTS},...
        'bayesFile2',tradespace.bayesFile2{idxTS},...
        'conflictTime_s',tradespace.conflictTime_s(idxTS,:),...
        'dofMaxRange_ft',tradespace.dofMaxRange_ft(idxTS,:),...
        'dofMaxVert_ft',tradespace.dofMaxVert_ft(idxTS,:),...
        'initHorz_ft',tradespace.initHorz_ft(idxTS,:),...
        'initVert_ft',tradespace.initVert_ft(idxTS,:),...
        'mdlType1',tradespace.mdlType1(idxTS,:),...
        'mdlType2',tradespace.mdlType2(idxTS,:),...
        'dem',dem,...
        'initialSeed',initialSeeds(ii),...
        'outHash',outHash,...
        'isGeoPointOverride',true,...
        'isSampleAlt1',true,...
        'isStartEndNonZero1',true,...
        'isPlot',true,...
        'isZip',false);
    
    % Record performance
    workTime_s = toc;
    nEncounters = size(metadata,1);
    
    % Save to status table
    status.workTime_s(ii) = workTime_s;
    status.nEncounters(ii) = nEncounters;
    
    save([outDir filesep 'metadata_' outHash '.mat'],'metadata','workTime_s','nEncounters','ii','iiAnchors');
end

%% Save Display status to screen
save([outDirRoot filesep 'tradespace.mat'],'status','-append');
warning('on','all')
disp('Done!');
