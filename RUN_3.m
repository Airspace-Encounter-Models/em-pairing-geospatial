% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
% Directory where encounters are stored
outDirRoot = [getenv('AEM_DIR_GEOPAIR') filesep 'output' filesep '2_encounters']

% Desired number of encounters
nEncWant = 1e6;

% Random seed
seed = 1;

%% Update directory where encounters are stored
if ispc
    outDirRoot = ['Z:' outDirRoot(2:end)];
end

%% Set random seed
rng(seed,'twister');

%% Preallocate
if isempty(bayesFile1)
    bfs = strings(numel(bayesFile2),2);
    bfs(:,1) = bayesFile1;
    bfs(:,2) = bayesFile2;
else
    bfs = allcomb(bayesFile1,bayesFile2);
end

tradespaceFilter = cell(size(bfs,1),1);
tasksFilter = cell(size(bfs,1),1);
listing = cell(size(bfs,1),1);
metadata = cell(size(bfs,1),1);
files = cell(size(bfs,1),1);

%% Load tradespace and tasks
load([outDirRoot filesep 'tradespace.mat'],'tradespace','tasks','anchors','initialSeeds','anchorRange_nm','usecase');

%% Iterate and identify
for ii=1:1:size(bfs,1)
    [tradespaceFilter{ii}, tasksFilter{ii}, listing{ii},metadata{ii}] = IdentifyEncounterFiles(tradespace,tasks,outDirRoot,...
        'bayesFile1',bfs{ii,1},'bayesFile2',bfs{ii,2},'nEncWant',nEncWant);
    files{ii} = string(compose('%s/enc_%s.dat',string(listing{ii}.folder(listing{ii}.isUse)),string(listing{ii}.nameHash(listing{ii}.isUse))));
    fprintf('%i/%i\n',ii,size(bfs,1));
end

%% Save
save([outDirRoot filesep 'listing.mat'],'tradespaceFilter','tasksFilter','listing','files','nEncWant','seed','bayesFile1','bayesFile2');
save([outDirRoot filesep 'metadata.mat'],'metadata','tasksFilter','listing','files','nEncWant','seed','bayesFile1','bayesFile2','-v7.3');
