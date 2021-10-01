function [tasks, anchors, seeds] = CreateTasksPairingGeo(tradespace,anchorRange_nm,usecase,seedGen)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Input handling
if nargin < 2; anchorRange_nm = 1.00; end
if nargin < 3; usecase = 'all'; end
if nargin < 4; seedGen = 0; end

%% Inputs hardcode
chunkSize = 5000;
maxSeed = 2^30; %

%% Set random seed
if ~isnan(seedGen) && ~isempty(seedGen)
    oldSeed = rng;
    rng(seedGen,'twister');
end

%% Prealocate
anchors = table.empty(0,6);

%% Load all anchors

switch usecase
    case 'uncor_uncor'
        
        inFile = [getenv('AEM_DIR_GEOPAIR') filesep 'output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat','uncor_uncor',anchorRange_nm)];
        
        if exist(inFile) == 2
             data = load(inFile,'anchors','anchorRange_nm');
             anchors = data.anchors;
        else
            % Generate if it doesn't exist
            iso_3166_2 = {'US-CO';'US-HI';'US-KS';'US-MA';'US-MS';'US-NC';'US-ND';'US-NV';'US-NY'};
            for ii=iso_3166_2'
                data = preallocAnchors( ii{1}, anchorRange_nm, '', false);
                data.num_geospatial = ones(size(data,1),1);
                data.iso_3166_2 = repmat(ii,size(data,1),1);
                
                if isempty(anchors);
                    anchors = data;
                else
                    anchors = [anchors;data];
                end
            end
            save(inFile,'anchors','iso_3166_2','anchorRange_nm');
        end
        
    otherwise
        uIso = string(unique(tradespace.subdivision1,'stable')');
        for ii=uIso
            % Create input filename containing the output of RUN_1
            inFile = [getenv('AEM_DIR_GEOPAIR') filesep 'output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat',sprintf('%s_%s',usecase,ii),anchorRange_nm)];
            
            % Load anchors
            data = load(inFile,'anchors','anchorRange_nm');
            data = data.anchors;
            data.iso_3166_2 = repmat(ii,size(data,1),1);
            
            if isempty(anchors);
                anchors = data(data.num_geospatial > 0,:);
            else
                anchors = [anchors; data(data.num_geospatial > 0,:)];
            end
        end 
end

% Sort by iso_3166_2 first (this is important for the next steps)
anchors = sortrows(anchors,{'iso_3166_2','num_geospatial'},{'ascend','descend'});
anchors.anchorRange_nm = repmat(anchorRange_nm,size(anchors,1),1);

%% Calculate start and end indicies for tasks w.r.t to anchors
nAnchors = size(anchors,1);
[uIso,~,ic] = unique(anchors.iso_3166_2,'stable');

agg.iso_3166_2 = strings(0,0);
agg.sidx = [];
agg.eidx = [];
for ii=1:1:numel(uIso)
    % Filter
    l = ic==ii;
    nii = nnz(l);
    B = cumsum(anchors.num_geospatial(l));
    
    % Calculate indicies
    idx = find(diff(mod(B,chunkSize)) <= 0);
    if isempty(idx); idx = 1; end
    if idx(1) ~= 1; idx = [1; idx]; end
    if idx(end) ~= nii; idx = [idx; nii]; end
    sidx = idx(1:end-1);
    eidx = idx(2:end)-1; eidx(end) = nii;
    
    % Assign
    agg.iso_3166_2 = [agg.iso_3166_2;repmat(uIso(ii),size(sidx))];
    if ii > 1
        agg.sidx = [agg.sidx; sidx + agg.eidx(end)];
        agg.eidx = [agg.eidx; eidx + agg.eidx(end)];
    else
        agg.sidx = sidx;
        agg.eidx = eidx;
    end
end

iso_3166_2 = agg.iso_3166_2;
%sidx = agg.sidx;
%eidx = agg.eidx;
clear tasks;

%% Concat start and end incides for tasks w.r.t to tradespace
configId = [];
sidx = [];
eidx = [];
for ii=1:1:size(tradespace,1)
    iiId = tradespace.configId(ii);
    
    switch usecase
        case 'uncor_uncor'
            isIso = true(size(agg.iso_3166_2));
        otherwise
            isIso = strcmpi(tradespace.subdivision1{ii},agg.iso_3166_2); 
    end
    nz = nnz(isIso);
    
    if ii > 1
        configId = [configId ; repmat(iiId,nz,1)];
        sidx = [sidx; agg.sidx(isIso)];
        eidx = [eidx; agg.eidx(isIso)];
    else
        configId = repmat(iiId,nz,1);
        sidx = agg.sidx(isIso);
        eidx = agg.eidx(isIso);
    end
end
globalId = (1:1:numel(configId))';

% Create task table and seed array
tasks = table(globalId,configId,sidx,eidx);
seeds = randi(maxSeed,size(tasks,1),1);

%% Change back to original seed
if ~isnan(seedGen) && ~isempty(seedGen)
    rng(oldSeed);
end

