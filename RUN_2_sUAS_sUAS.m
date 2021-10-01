% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
tradespace = CreateTradespacePairingGeo('mdlType1','geospatial','mdlType2','geospatial',...
    'initHorz_ft',[4405 8810],'initVert_ft',[0 300],...
    'thresHorz_ft',0.6 * unitsratio('ft','nm'),'thresVert_ft',100); %https://arxiv.org/abs/1911.00110

% Parent output directory
outDirRoot = [getenv('AEM_DIR_GEOPAIR') filesep 'output' filesep '2_encounters' filesep 'geospatial_uncor' filesep datestr(date,'yyyy-mm-dd')];

% RUN_1 parameters
anchorRange_nm = 1.00;
usecase = 'all';

%% Suppress warning
warning('off','map:io:UnableToDetermineCoordinateSystemType')

%% Save tradespace
[~, ~, ~] = mkdir(outDirRoot);
save([outDirRoot filesep 'tradespace.mat'],'tradespace');

%% Create status table to record performance
status = table( tradespace.configId,zeros(size(tradespace,1),1),zeros(size(tradespace,1),1),'VariableNames',{'tradespaceId','runTime_s','nEncounters'});

%% Iterate over tradespace
for ii=1:1:size(tradespace,1)
    tic
    % Create input filename containing the output of RUN_1
    inFile = ['output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat',sprintf('%s_%s',usecase,tradespace.subdivision1{ii}),anchorRange_nm)];
    
    % Create output directory
    outDir = [outDirRoot filesep num2str(ii)];
    [~, ~, ~] = mkdir(outDir);
    
    % Create encounters
    metadata = createEncounters_2(inFile,outDir,...
        tradespace.encTime_s(ii,:),tradespace.thresHorz_ft(ii,:),tradespace.thresVert_ft(ii,:),...
        'bayesFile1',tradespace.bayesFile1{ii},...
        'bayesFile2',tradespace.bayesFile2{ii},...
        'conflictTime_s',tradespace.conflictTime_s(ii,:),...
        'dofMaxRange_ft',tradespace.dofMaxRange_ft(ii,:),...
        'dofMaxVert_ft',tradespace.dofMaxVert_ft(ii,:),...
        'initHorz_ft',tradespace.initHorz_ft(ii,:),...
        'initVert_ft',tradespace.initVert_ft(ii,:),...
        'mdlType1',tradespace.mdlType1(ii,:),...
        'mdlType2',tradespace.mdlType2(ii,:),...
        'isStartEndNonZero1',true,...
        'isStartEndNonZero2',true,...
        'isPlot',true);
    
    status.runTime_s(ii) = toc;
    status.nEncounters(ii) = size(metadata,1);
    save([outDirRoot filesep sprintf('metadata_%05.f',ii)],'metadata','status','inFile','ii','tradespace');
end

%% Save Display status to screen
save([outDirRoot filesep 'tradespace.mat'],'status','-append');
warning('on','all')
disp('Done!');
