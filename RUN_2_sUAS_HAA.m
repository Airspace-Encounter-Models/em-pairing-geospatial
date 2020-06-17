%% Inputs
iso_3166_2 = {'US-CA','US-FL','US-KS','US-MA','US-MS','US-NC','US-ND','US-NH','US-NV','US-NY','US-OK','US-PR','US-RI','US-TN','US-TX','US-VA'}; 
%iso_3166_2 = {'US-MA','US-MS','US-NH','US-NY','US-PR','US-TN'}; % Smaller states
%iso_3166_2 = {'US-VA','US-CA','US-FL','US-KS','US-ND','US-NV','US-OK','US-TX','US-NC'}; % Larger states
%iso_3166_2 = {'US-RI'}; % Smallest state

% RUN_1 parameters
anchorRange_nm = 1.25;
usecase = 'shield';

% Encounters
encTime_s = 60; % Total encounter time, where CPA occurs at encTime_s /2
initHorz_ft = [6076 round(unitsratio('ft','nm') * anchorRange_nm * 2)]; % [1 2.5] nautical miles
initVert_ft = [0 500];
thresHorz_ft = 500; % CPA <= minHorz_ft % Conflict threshold criteria
thresVert_ft = 100;

% Anchors
anchorPercent = 1;
maxEncPerPair = 20;

% Bayes
dofMaxRange_ft = 200;
z_agl_tol_ft = 200;
demBayes = 'dted1';

% Altitude Sampling
if strcmp(usecase,'shield')
    isSampleAlt1 = false;
    maxAlt1_ft_agl = [];
else
    isSampleAlt1 = true;
    maxAlt1_ft_agl = 1200;
end
endisSampleAlt2 = false;

% Airspeed Sampling
vrange1_kts = [1 10];

%% HAA trajectories
inDir2 = arrayfun(@(x)([getenv('AEM_DIR_HAA') filesep 'output' filesep 'v1-2020-01-14-HAA' filesep num2str(x) 'ft']),0:1e2:500,'uni',false)';

%% Create status table to record performance
if ~iscolumn(iso_3166_2)
    iso_3166_2 = iso_3166_2';
end
status = table(iso_3166_2,zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),'VariableNames',{'iso_3166_2','runtime_s','n_encounters','n_pairs'});
outHash = ['_' usecase '_haa' '_anchorPercent' num2str(anchorPercent * 100) '_maxEncPerPair' num2str(maxEncPerPair) '_encTime' num2str(encTime_s) '_thresHorz_ft' num2str(thresHorz_ft) '_thresVert_ft' num2str(thresVert_ft) '_inithorzmin' num2str(initHorz_ft(1)) '_inithorzmax' num2str(initHorz_ft(2)) '_initvertmin' num2str(initVert_ft(1)) '_initvertmax' num2str(initVert_ft(2))];

%% Iterate through adminstrative boundaries
for i=1:1:numel(iso_3166_2)
    tic
    % Create input filename
    inFile = ['output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat',sprintf('%s_%s',usecase,iso_3166_2{i}),anchorRange_nm)];
    
    % Create output directory
    outDir = ['output' filesep 'NMAC' filesep iso_3166_2{i} outHash];
    
    % Create encounters
    [status.n_encounters(i), anchors, files1, files2, isEncounter]  = createEncounters_2(inFile,encTime_s,thresHorz_ft,thresVert_ft,outDir,...
        'acmodel1','geospatial',...
        'acmodel2','bayes',...
        'inDir2',inDir2,...
        'maxEncPerPair',10,...
        'initHorz_ft',initHorz_ft,...
        'initVert_ft',initVert_ft,...
        'isSampleAlt1',isSampleAlt1,...
        'isSampleAlt2',isSampleAlt2,...
        'vrange1_kts',vrange1_kts,...
        'maxAlt1_ft_agl',maxAlt1_ft_agl,...
        'demBayes',demBayes,...
        'dofMaxRange_ft',dofMaxRange_ft,...
        'z_agl_tol_ft',z_agl_tol_ft,...
        'rowstart2',1,...
        'anchorPercent',anchorPercent,...
        'maxEncPerPair',maxEncPerPair,...
        'labelZ','barometicAlt_ft_msl',...
        'isPlot',true);
    
    status.runtime_s(i) = round(toc);
    status.n_pairs(i) = sum(isEncounter);
    
    save(['output' filesep 'NMAC' filesep '2_encounters_' iso_3166_2{i} outHash],'anchors','files1','files2','isEncounter','thresHorz_ft','encTime_s','anchorPercent','maxEncPerPair');
end

%% Save Display status to screen
save(['output' filesep 'NMAC' filesep '2_performanace' outHash],'status','anchorRange_nm','usecase','encTime_s','initHorz_ft','initVert_ft','thresHorz_ft','isSampleAlt1','isSampleAlt2','anchorPercent','maxEncPerPair');
disp('Done!');