function tradespace = CreateTradespacePairingGeo(varargin)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO createEncounters_2

%% Input parser
if nargin < 1; varargin = cell.empty(0,1); end
p = inputParser;

% Initial Seperation
addParameter(p,'initHorz_ft',[6076 97217],@isnumeric); %97217 ~ 16 nm ~ 421 feet (250 knots) * 240 seconds (uncor MVL)
addParameter(p,'initVert_ft',[0 1500],@isnumeric);

% Conflict Threshold
addParameter(p,'thresHorz_ft',[500; 2000; 4500],@iscolumn);
addParameter(p,'thresVert_ft',[100; 250; 450],@iscolumn);

% Time
addParameter(p,'encTime_s',[90; 150; 210],@iscolumn);
addParameter(p,'conflictTime_s',[60; 120; 180],@iscolumn);

% Obstacles
addParameter(p,'dofMaxRange_ft',[500],@isnumeric);
addParameter(p,'dofMaxVert_ft',[500],@isnumeric); % Minimum height of obstacles

% Aircraft types
addParameter(p,'mdlType1','geospatial',@(x) ischar(x) && any(strcmpi(x,{'bayes','geospatial'})));
addParameter(p,'mdlType2','bayes',@(x) ischar(x) && any(strcmpi(x,{'bayes','geospatial'})));
defaultBayes = compose('%s/%s',[getenv('AEM_DIR_BAYES') filesep 'model'],["uncor_1200exclude_rotorcraft_v1p2.txt"; "uncor_1200only_rotorcraft_v1p2.txt";"uncor_1200exclude_fwse_v1p2.txt"; "uncor_1200only_fwse_v1p2.txt";"uncor_1200exclude_fwme_v1p2.txt"; "uncor_1200only_fwme_v1p2.txt"]);
addParameter(p,'bayesFile1',defaultBayes);
addParameter(p,'bayesFile2',defaultBayes);
addParameter(p,'iso_3166_2',{'US-CO';'US-HI';'US-KS';'US-MA';'US-MS';'US-NC';'US-ND';'US-NV';'US-NY'});

% Parse
parse(p,varargin{:});
parms = p.Results;
fields = fieldnames(parms);

%% Input Handling
assert(all(p.Results.encTime_s >= p.Results.conflictTime_s),'Failed to satisfy encTime_s >= conflictTime_s');

%% Initial Seperation, Conflict, and Time Thresholds
% Number of elements for each variable / group
nInit = size(parms.initHorz_ft,1); % Initial seperation
nThres = size(parms.thresHorz_ft,1); % Conflict seperation threshold
nTime = size(parms.encTime_s,1); % Encounter and conflict time s
nObs = size(parms.dofMaxRange_ft,1); % Obstacles

% Combination of variables
A = allcomb(1:1:nInit,1:1:nThres,1:1:nTime,1:1:nObs,'matlab');

% Create tradespace with initial and threshold combinations
tradespace = table(parms.initHorz_ft(A(:,1),:),parms.initVert_ft(A(:,1),:),...
    parms.thresHorz_ft(A(:,2),:), parms.thresVert_ft(A(:,2),:),...
    parms.encTime_s(A(:,3),:), parms.conflictTime_s(A(:,3),:),...
    parms.dofMaxRange_ft(A(:,4),:),parms.dofMaxVert_ft(A(:,4),:),...
    'VariableNames',{'initHorz_ft','initVert_ft','thresHorz_ft','thresVert_ft','encTime_s','conflictTime_s','dofMaxRange_ft','dofMaxVert_ft'});

%% Aircraft Types
switch parms.mdlType1
    case 'bayes'
        bayes1 = parms.bayesFile1;
        geo1 = char(32); % Unicode whitespace
    case 'geospatial'
        geo1 = parms.iso_3166_2;
        bayes1 = string;
end

switch parms.mdlType2
    case 'bayes'
        bayes2 = parms.bayesFile2;
        geo2 = char(32); % Unicode whitespace
    case 'geospatial'
        geo2 = parms.iso_3166_2;
        bayes2 = string;
end

nb1 = size(bayes1,1); nb2 = size(bayes2,1);
ng1 = size(geo1,1); ng2 = size(geo2,1);

% Combination of variables
% bayes1, bayes2, geo1, geo2
A = allcomb(1:1:nb1,1:1:nb2,1:1:ng1,1:1:ng2,'matlab');

% Iterate over combinations
oldTradespace = tradespace;
for ii=1:1:size(A,1)
    % Get unappended table
    iiT = oldTradespace;
    nt = size(iiT,1);
    
    % Append table
    iiT.mdlType1 = repmat(parms.mdlType1,nt,1);
    iiT.mdlType2 = repmat(parms.mdlType2,nt,1);
    
    iiT.bayesFile1 = repmat(bayes1(A(ii,1),:),nt,1);
    iiT.bayesFile2 = repmat(bayes2(A(ii,2),:),nt,1);
    iiT.subdivision1 = repmat(geo1(A(ii,3),:),nt,1);
    iiT.subdivision2 = repmat(geo2(A(ii,4),:),nt,1);
    
    % Concat
    if ii == 1;
        tradespace = iiT;
    else
        tradespace = [tradespace ; iiT];
    end
end

%% Numerical id
tradespace.configId = (1:1:size(tradespace,1))';
tradespace = movevars(tradespace,'configId','Before',1);

