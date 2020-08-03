% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
%iso_3166_2 = {'US-MA'};%,'US-CA','US-TX'};%'US-FL','US-KS','US-MA','US-MS','US-NC','US-ND','US-NH','US-NV','US-NY','US-OK','US-PR','US-RI','US-TN','US-TX','US-VA'};
iso_3166_2 = {'US-KS','US-MA','US-MS','US-NC','US-ND','US-NV','US-TN','US-VA'};%,'US-CA','US-TX'};%'US-FL','US-KS','US-MA','US-MS','US-NC','US-ND','US-NH','US-NV','US-NY','US-OK','US-PR','US-RI','US-TN','US-TX','US-VA'};
%iso_3166_2 = {'US-MA','US-MS','US-NH','US-NY','US-PR','US-TN'}; % Smaller states
%iso_3166_2 = {'US-VA','US-CA','US-FL','US-KS','US-ND','US-NV','US-OK','US-TX','US-NC'}; % Larger states
%iso_3166_2 = {'US-RI'}; % Smallest state

% RUN_1 parameters
anchorRange_nm = 1.00;
usecase = 'all-noshield';

% Encounters
encTime_s = 60; % Total encounter time, where CPA occurs at encTime_s /2
initHorz_ft = [6076 round(unitsratio('ft','nm') * anchorRange_nm * 2)]; % [1 2.5] nautical miles
initVert_ft = [0 500];
thresHorz_ft = 500; % CPA <= minHorz_ft % Conflict threshold criteria
thresVert_ft = 100;

% Anchors
anchorPercent = 1;
maxEncPerPair = 2;

% Altitude and Airspeed Sampling
if strcmp(usecase,'shield')
    isSampleAlt1 = false;
    maxAlt1_ft_agl = []; % Not used because isSampleAlt1 = false
    vrange1_kts = [1 10];
else
    isSampleAlt1 = true;
    maxAlt1_ft_agl = 700;
    vrange1_kts = [10 87];
end
minAlt1_ft_agl = 100;
isSampleAlt2 = false;

% Bayes
bayesModel = 'uncor_allcode_rotorcraft_v1'; % Model name
bayesDate = '31-Jul-2020'; % Date created
minAlt2_ft_agl = minAlt1_ft_agl;
maxAlt2_ft_agl = round(700 + max(initVert_ft),-2); % The ASTM sUAS TOR has a ceiling of 1200 feet AGL
num_files2 = 40; % Number of files dependent upon em-model-manned-bayes/RUN_1_emsample
dofMaxRange_ft = 250; % Maximum "size" of obstacles, a circle around the obstacle will be generated with radius = dofMaxRange_ft
z_agl_tol_ft = 200;
demBayes = 'dted1';

% Airspace class char
% 1 = B = Class B
% 2 = C = Class C
% 3 = D = Class D
% 4 = O = Other (Class A, E, G)
classInclude = {'O'};

%% Create status table to record performance
if ~iscolumn(iso_3166_2)
    iso_3166_2 = iso_3166_2';
end
status = table(iso_3166_2,zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),'VariableNames',{'iso_3166_2','runtime_s','n_encounters','n_pairs'});

%% Iterate through adminstrative boundaries
for i=1:1:numel(iso_3166_2)
    tic
    % Create input filename
    inFile = ['output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat',sprintf('%s_%s',usecase,iso_3166_2{i}),anchorRange_nm)];
    
    % Set random seed
    rng(i);
    
    % Iterate over airspace classes
    for j=classInclude
        % Geographic string
        % Puerto Rico and Hawaii are islands
        if any(strcmp(iso_3166_2{i},{'US-PR','US-HI'}))
            G = 'G2'; % Islands
        else
            G = 'G1'; % CONUS + AK
        end
        
         % Airspace Class String
        switch j{1}
            case 'B'
                A = ['A' '1'];
            case 'C'
                A = ['A' '2'];
            case 'D'
                A = ['A' '3'];
            case 'O'
                A = ['A' '4'];
            otherwise
                error('classInclude not recognized\n');
        end
             
        % Uncorrelated trajectories
        % Subdirectories are organized by geographic variable (G), airspace class (A), and altitude layer (L)
        inDirG = arrayfun(@(x)([getenv('AEM_DIR_BAYES') filesep 'output' filesep 'tracks' filesep bayesModel '_' bayesDate filesep num2str(x)]),randperm(num_files2,round(num_files2/2)),'uni',false); % randperm(num_files2,20) is used to reduce the number of files searched for in createEncounters_2
        inDirGA = cellfun(@(x)([x filesep G filesep A]),inDirG,'uni',false);
        inDir2 = {};
        for k=1:1:numel(inDirGA)
            inDirK = arrayfun(@(x)([inDirGA{k} filesep num2str(x) 'ft']),minAlt2_ft_agl:1e2:maxAlt2_ft_agl,'uni',false)'; % The altitude range corresponds to the directory structure from em-model-manned-bayes/RUN_2_sample2track
            inDir2 = [inDir2; inDirK];
        end
        
        % Create output hash
        outHash = ['_' usecase '_' strrep(bayesModel,'_','-') '_' G '_' A '_anchorPercent' num2str(anchorPercent * 100) '_maxEncPerPair' num2str(maxEncPerPair) '_encTime' num2str(encTime_s) '_thresHorz_ft' num2str(thresHorz_ft) '_thresVert_ft' num2str(thresVert_ft) '_inithorzmin' num2str(initHorz_ft(1)) '_inithorzmax' num2str(initHorz_ft(2)) '_initvertmin' num2str(initVert_ft(1)) '_initvertmax' num2str(initVert_ft(2))];
        
        % Create output directory
        outDir = ['output' filesep 'NMAC' filesep 'AC1speed' num2str(vrange1_kts(1)) '-' num2str(vrange1_kts(2)) filesep  iso_3166_2{i} outHash];
        
        % Create encounters
        [status.n_encounters(i), anchors, files1, files2, isEncounter] = createEncounters_2(inFile,encTime_s,thresHorz_ft,thresVert_ft,outDir,...
            'acmodel1','geospatial',...
            'acmodel2','bayes',...
            'inDir2',inDir2,...
            'maxEncPerPair',10,...
            'initHorz_ft',initHorz_ft,...
            'initVert_ft',initVert_ft,...
            'isSampleAlt1',isSampleAlt1,...
            'isSampleAlt2',isSampleAlt2,...
            'vrange1_kts',vrange1_kts,...
            'minAlt1_ft_agl',minAlt1_ft_agl,...
            'maxAlt1_ft_agl',maxAlt1_ft_agl,...
            'demBayes',demBayes,...
            'dofMaxRange_ft',dofMaxRange_ft,...
            'z_agl_tol_ft',z_agl_tol_ft,...
            'rowstart2',1,...
            'anchorPercent',anchorPercent,...
            'maxEncPerPair',maxEncPerPair,...
            'classInclude',j,...
            'isPlot',false);
        
        save(['output' filesep 'NMAC' filesep '2_encounters_' iso_3166_2{i} outHash],'anchors','files1','files2','isEncounter','thresHorz_ft','encTime_s','anchorPercent','maxEncPerPair');       
    end
    status.runtime_s(i) = round(toc);
    status.n_pairs(i) = sum(isEncounter);
end

%% Save Display status to screen
save(['output' filesep 'NMAC' filesep '2_performanace' outHash],'status','anchorRange_nm','usecase','encTime_s','initHorz_ft','initVert_ft','thresHorz_ft','isSampleAlt1','isSampleAlt2','anchorPercent','maxEncPerPair','bayesDate','bayesModel');
disp('Done!');