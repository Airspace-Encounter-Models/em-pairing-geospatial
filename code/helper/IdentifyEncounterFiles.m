function [tsOut, tasksOut, listing, metadata] = IdentifyEncounterFiles(tsIn,tasksIn,inDir,varargin)

%% Input parser
if nargin < 1; varargin = cell.empty(0,1); end
p = inputParser;

% Required
addRequired(p,'tsIn',@istable);
addRequired(p,'tasksIn',@istable);
addRequired(p,'inDir',@ischar);

% Initial Seperation
%addParameter(p,'initHorz_ft',[6076 97217],@isnumeric); %97217 ~ 16 nm ~ 421 feet (250 knots) * 240 seconds (uncor MVL)
%addParameter(p,'initVert_ft',[0 1000],@isnumeric);

% Conflict Threshold
addParameter(p,'thresHorz_ft',[],@isnumeric);
addParameter(p,'thresVert_ft',[],@isnumeric);

% Time
addParameter(p,'encTime_s',[],@isnumeric);
addParameter(p,'conflictTime_s',[],@isnumeric);

% Obstacles
addParameter(p,'dofMaxRange_ft',[],@isnumeric);
addParameter(p,'dofMaxVert_ft',[],@isnumeric); % Minimum height of obstacles

% Aircraft types
addParameter(p,'bayesFile1',char.empty(0,0));
addParameter(p,'bayesFile2','uncor_1200only_fwse_v1p2');

% Filtering
addParameter(p,'nEncWant',1e6,@isnumeric);

% Parse
parse(p,tsIn,tasksIn,inDir,varargin{:});
parms = p.Results;

%%
tsOut = tsIn;
[~,tsOut.bayesFile1,~] = fileparts(tsIn.bayesFile1);
[~,tsOut.bayesFile2,~] = fileparts(tsIn.bayesFile2);

%% Logical filtering of tradespace

parNames = fields(parms);
varNames = tsOut.Properties.VariableNames;

tf = true(size(tsOut,1),numel(parNames));

for ii=1:1:numel(parNames)
    x = parms.(parNames{ii});
    if any(strcmp(parNames{ii},varNames)) && ~isempty(x)
        switch class(x)
            case 'char'
                tf(:,ii) = string(tsOut.(parNames{ii})) == string(x);
            otherwise
                tf(:,ii) = tsOut.(parNames{ii}) == x;
        end
    end
end

tf = all(tf,2);

tsOut = tsOut(tf,:);
tasksOut = tasksIn(ismember(tasksIn.configId,tsOut.configId),:);

%%
configId = [];
for ii=1:1:size(tsOut,1)
    iiId = tsOut.configId(ii);
    iiDir = sprintf('%s/%04.f',inDir,iiId);
    
    iiList = dir([iiDir filesep 'metadata*']);
    n = numel(iiList);
    
    if ii == 1
        listing = iiList;
        configId = repmat(iiId,n,1);
    else
        configId = [configId ; repmat(iiId,n,1)];
        listing = [listing; iiList];
    end
end

% Convert to table
listing = struct2table(listing);

% Parse name
listing.nameHash = strrep(strrep(listing.name,'metadata_',''),'.mat','');
listing.configId = configId;
C = cellfun(@(x)(strsplit(x,'_')),listing.nameHash,'UniformOutput',false);
listing.globalId = cellfun(@(x)(str2double(x{2})),C,'UniformOutput',true);

% Filter tasks
lg = ismember(tasksOut.globalId,listing.globalId);
lc = ismember(tasksOut.configId,listing.configId);
tasksOut = tasksOut(lg & lc,:);

%% Number of encounters and workTime_s
nEncounters = zeros(size(listing,1),1);
workTime_s = zeros(size(listing,1),1);
metadata = table;

for ii=1:1:size(nEncounters,1)
    % Load
    x = load([listing.folder{ii} filesep listing.name{ii}],'nEncounters','workTime_s','metadata');
    
    % Assign
    nEncounters(ii) = x.nEncounters;
    workTime_s(ii) = x.workTime_s;
    
    if ii==1;
        metadata = x.metadata;
    else
        metadata = [metadata; x.metadata];
    end
    
end
listing.nEncounters = nEncounters;
listing.workTime_s = workTime_s;

% tradespace
tsOut.nEncounters = nan(size(tsOut,1),1);
tsOut.workTime_s = nan(size(tsOut,1),1);

[uCId,~,~] = unique(listing.configId);
tsOut.nEncounters(ismember(tsOut.configId,listing.configId)) = arrayfun(@(x)(sum(listing.nEncounters(listing.configId == x))),uCId);
tsOut.workTime_s(ismember(tsOut.configId,listing.configId)) = arrayfun(@(x)(sum(listing.workTime_s(listing.configId == x))),uCId);

%% Filtering to desired number of encounters
% Random permutation of indicies
idx = randperm(size(listing,1));

% Preallocate
isUse = false(size(listing,1),1);
cIdx = 1;
cEnc = 0;

while cEnc < parms.nEncWant && cIdx < numel(idx)
    isUse(idx(cIdx)) = true;
    cEnc = cEnc + listing.nEncounters(idx(cIdx));
    cIdx = cIdx + 1;
end

listing.isUse = isUse;

