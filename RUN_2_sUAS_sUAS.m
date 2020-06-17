%% Inputs
iso_3166_2 = {'US-CA','US-FL','US-KS','US-MA','US-MS','US-NC','US-ND','US-NH','US-NV','US-NY','US-OK','US-PR','US-RI','US-TN','US-TX','US-VA'}; % dev cases
%iso_3166_2 = {'US-RI'}; % dev cases
anchorRange_nm = 0.17; % Distance pairs need to be away from anchor point (0.17 = unitsratio('nm','ft') * 1000)
pairsNum = 10000;
usecase = 'longlinearinfrastructure';

% Encounters
encTime_s = 60; % Total encounter time, where CPA occurs at encTime_s /2
initHorz_ft = [4405 8810];
initVert_ft = [0 300];
timeStep_s = 1;

% CPA threshold criteria
thresHorz_ft = 0.6 * unitsratio('ft','nm'); %https://arxiv.org/abs/1911.00110

%% Iterate through adminstrative boundaries
encCount = zeros(size(iso_3166_2));
for i=1:1:numel(iso_3166_2)
    % Create input filename
    inFile = ['output' filesep sprintf('pairs-%s-%0.2f-%i.mat',sprintf('%s-%s',usecase,iso_3166_2{i}),anchorRange_nm,pairsNum)];
    
    % Create output directory
    outDir = ['output' filesep iso_3166_2{i} '-' usecase '-' usecase '-cpa' num2str(thresHorz_ft) '-encTime' num2str(encTime_s) '-timestep' num2str(timeStep_s) '-inithorzmin' num2str(initHorz_ft(1)) '-inithorzmax' num2str(initHorz_ft(2)) '-initvertmin' num2str(initVert_ft(1)) '-initvertmax' num2str(initVert_ft(2))];
    
    % Create encounters
    encCount(i) = createEncounters_2(inFile,encTime_s,thresHorz_ft,outDir,...
        'maxEncPerPair',10,...
        'acmodel1','geospatial',...
        'acmodel2','geospatial',...
        'initHorz_ft',initHorz_ft,...
        'initVert_ft',initVert_ft,...
        'timeStep_s', timeStep_s,...
        'isPlot',false);
end

%% Display status to screen
disp('Done!');