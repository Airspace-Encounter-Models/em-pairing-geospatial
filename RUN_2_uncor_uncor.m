% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
% Tradespace
tradespace = CreateTradespacePairingGeo('mdlType1','bayes','mdlType2','bayes',...
    'initHorz_ft',[12152 97217],'initVert_ft',[0 1500],'thresHorz_ft',[4000],'thresVert_ft',450,...
    'conflictTime_s',[60; 120; 180],'encTime_s',[90; 150; 210],...
    'dofMaxRange_ft',500,'dofMaxVert_ft',500);

% Parent output directory
outDirRoot = ['output' filesep '2_encounters' filesep 'uncor_uncor' filesep datestr(date,'yyyy-mm-dd')];

% Mock RUN_1_ parameters
anchorRange_nm = 1.00;

% Other parameter not defined in tradespace table
rangeAlt1_ft_agl = [50 1200];
rangeAlt2_ft_agl = [50 2000];
dem = 'globe';
maxTries = 12;
seedGen = 1;

%% Suppress warning
warning('off','map:io:UnableToDetermineCoordinateSystemType')

%% Create tasks
[tasks, anchors,initialSeeds] = CreateTasksPairingGeo(tradespace,anchorRange_nm,'uncor_uncor',seedGen);

%% Save tradespace
[~, ~, ~] = mkdir(outDirRoot);
save([outDirRoot filesep 'tradespace.mat'],'tradespace');

%% Create status table to record performance
status = table( tradespace.configId,zeros(size(tradespace,1),1),zeros(size(tradespace,1),1),'VariableNames',{'tradespaceId','runTime_s','nEncounters'});

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
        'maxTries',maxTries,...
        'rangeAlt1_ft_agl',rangeAlt1_ft_agl,...
        'rangeAlt2_ft_agl',rangeAlt2_ft_agl,...
        'isSampleAlt1',false,...
        'isSampleAlt2',false,...
        'isSampleV1',false,...
        'isSampleV2',false,...
        'isStartEndNonZero1',false,...
        'isStartEndNonZero2',false,...
        'isPlot',false,...
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
