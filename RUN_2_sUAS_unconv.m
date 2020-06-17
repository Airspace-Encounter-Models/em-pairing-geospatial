%% Inputs
iso_3166_2 = {'US-CA','US-FL','US-KS','US-MA','US-MS','US-NC','US-ND','US-NH','US-NY','US-NV','US-OK','US-PR','US-RI','US-TN','US-TX','US-VA'}; % dev cases
%iso_3166_2 = {'US-MA','US-MS','US-NH','US-NY','US-PR','US-TN'}; % Smaller states
%iso_3166_2 = {'US-VA','US-CA','US-FL','US-KS','US-ND','US-NV','US-OK','US-TX','US-NC'}; % Larger states
%iso_3166_2 = {'US-RI'}; % Smallest state

% RUN_1 parameters
anchorRange_nm = 1.25;
usecase = 'all-noshield';

% Encounters
encTime_s = 60; % Total encounter time, where CPA occurs at encTime_s /2
initVert_ft = [0 500];
thresHorz_ft = 500; % CPA <= minHorz_ft % Conflict threshold criteria
thresVert_ft = 100;

% Anchors
anchorPercent = 1;
maxEncPerPair = 20;

% Bayes
z_agl_tol_ft = 200;
demBayes = 'dted1';
models = {'glider_v1p2_08-Mar-2020','paraglider_v1p2_08-Mar-2020'};
model = models{1}; % Manually change

switch model
    case 'glider_v1p2_08-Mar-2020'
        dofMaxRange_ft = 500;
        initHorz_ft = [6076 round(unitsratio('ft','nm') * anchorRange_nm * 2)]; % [1 2.5] nautical miles
        maxHeightBayes_ft = 1700;
        aclass = {'O'};
    case 'paraglider_v1p2_08-Mar-2020'
        dofMaxRange_ft = 200;
        initHorz_ft = [1000 round(unitsratio('ft','nm') * anchorRange_nm * 2)]; % [1 2.5] nautical miles
        maxHeightBayes_ft = 1500;
        aclass = {'B','C','D','O'};
end

% Altitude and Airspeed Sampling
if strcmp(usecase,'shield')
    isSampleAlt1 = false;
    maxAlt1_ft_agl = [];
    vrange1_kts = [1 10];
else
    isSampleAlt1 = true;
    maxAlt1_ft_agl = 1200;
    vrange1_kts = [10 60];
end
isSampleAlt2 = false;

%% Manned trajectories
dirBayes = [getenv('AEM_DIR_BAYES') filesep 'output' filesep 'tracks' filesep model];
listing = dir(dirBayes);
listing(ismember( {listing.name}, {'.', '..'})) = [];  %remove . and ..

%% Create status table to record performance
if ~iscolumn(iso_3166_2)
    iso_3166_2 = iso_3166_2';
end
status = table(iso_3166_2,zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),'VariableNames',{'iso_3166_2','runtime_s','n_encounters','n_pairs'});
outHash = ['_' usecase '_unconv_' model '_anchorPercent' num2str(anchorPercent * 100) '_maxEncPerPair' num2str(maxEncPerPair) '_encTime' num2str(encTime_s) '_thresHorz_ft' num2str(thresHorz_ft) '_thresVert_ft' num2str(thresVert_ft) '_inithorzmin' num2str(initHorz_ft(1)) '_inithorzmax' num2str(initHorz_ft(2)) '_initvertmin' num2str(initVert_ft(1)) '_initvertmax' num2str(initVert_ft(2))];

%% Iterate through adminstrative boundaries
for i=1:1:numel(iso_3166_2)
    tic
    
    % Set random seed
    rng(i);
    
    % Create input filename
    inFile = ['output' filesep sprintf('pairs-%s-anchorRange%0.2f.mat',sprintf('%s_%s',usecase,iso_3166_2{i}),anchorRange_nm)];
    
    % Create output directory
    outDir = ['output' filesep 'NMAC' filesep 'AC1speed' num2str(vrange1_kts(1)) '-' num2str(vrange1_kts(2)) filesep  iso_3166_2{i} outHash];
    
    % Bayes
    inDirBayes = arrayfun(@(x)([dirBayes filesep num2str(x)]),randperm(numel(listing),10),'uni',false)';
    inDir2 = {};
    for k=1:1:numel(inDirBayes)
        inDirK = arrayfun(@(x)([inDirBayes{k} filesep num2str(x) 'ft']),0:1e2:maxHeightBayes_ft,'uni',false)';
        inDir2 = [inDir2; inDirK];
    end
    
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
        'maxAlt1_ft_agl',maxAlt1_ft_agl,...
        'demBayes',demBayes,...
        'dofMaxRange_ft',dofMaxRange_ft,...
        'z_agl_tol_ft',z_agl_tol_ft,...
        'rowstart2',1,...
        'vrange1_kts',vrange1_kts,...
        'anchorPercent',anchorPercent,...
        'maxEncPerPair',maxEncPerPair,...
        'classInclude',aclass,...
        'isPlot',false);
    
    status.runtime_s(i) = round(toc);
    status.n_pairs(i) = sum(isEncounter);
    
    save(['output' filesep 'NMAC' filesep '2_encounters_' iso_3166_2{i} outHash],'anchors','files1','files2','isEncounter','thresHorz_ft','encTime_s','anchorPercent','maxEncPerPair');
end

%% Save Display status to screen
save(['output' filesep 'NMAC' filesep '2_performanace' outHash],'status','anchorRange_nm','usecase','encTime_s','initHorz_ft','initVert_ft','thresHorz_ft','isSampleAlt1','isSampleAlt2','anchorPercent','maxEncPerPair');
disp('Done!');