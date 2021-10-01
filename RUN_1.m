% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% INPUTS
% Adminstrative Boundaries and Trajectory Parameters
iso_3166_2 = {'US-CO';'US-HI';'US-KS';'US-MA';'US-MS';'US-NC';'US-ND';'US-NV';'US-NY'};

% Pairs parameters
% 30*(150 / 3600) = 1.25...closing speed of 150 knots for 30 seconds
% 2000 ft = well clear
anchorRange_nm = 1.00; % Distance pairs need to be away from anchor point

%% Inputs hardcode
usecase = 'all';

% Root directory of geospatial trajectories
inDirRoot = ['~', filesep 'CASSATT_shared' filesep 'geospatial_trajectories'];

%% Create status table to record performance
if ~iscolumn(iso_3166_2)
    iso_3166_2 = iso_3166_2';
end
status = table(iso_3166_2,zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),zeros(size(iso_3166_2)),'VariableNames',{'iso_3166_2','runtime_s','nPoints','nGeoFiles','nClusters'});

%% Iterate through adminstrative boundaries
for ii=1:1:numel(iso_3166_2)
    tic
    
    % Parse directories of the geospatial trajectories
    [inDir,~] = parseGeoTrajDirectory([inDirRoot filesep iso_3166_2{ii}], usecase);
    
    if isempty(inDir)
        fprintf('No input directories for iso_3166_2 = %s\n',iso_3166_2{ii});
    else
        % Execute
        [anchors, inFiles, numClust] = findPairs_1(iso_3166_2{ii},inDir,sprintf('%s_%s',usecase,iso_3166_2{ii}),'anchorRange_nm',anchorRange_nm,'seed',ii);
        
        % Record stats
        status.runtime_s(ii) = round(toc);
        status.nPoints(ii) = size(anchors,1);
        status.nGeoFiles(ii) = numel(inFiles);
        status.nClusters(ii) = numClust;
        
        % These are saved and we only outputted for performance, don't need
        % them hanging around for the next loop
        clear anchors inFiles
    end
end

%% Save Display status to screen
save(['output' filesep '1_performanace_' usecase '_' 'anchorRange' num2str(anchorRange_nm) '_' date '.mat'],'status','anchorRange_nm','usecase');
disp('Done!');

