function [output,listing] = parseGeoTrajDirectory(inDir, usecase)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

% Input handling
if nargin < 1; inDir = [getenv('AEM_DIR_GEOSPATIAL') filesep 'output' filesep 'trajectories']; end
if nargin < 2; usecase = 'all'; end

% Input directories of geospatial trajectories
listing = dir([inDir '*' filesep '*']);
listing(ismember( {listing.name}, {'.', '..'})) = [];  %remove . and ..
isShield = contains({listing.folder},'shield');
output = cellfun(@(f,n)([f filesep n]),{listing.folder},{listing.name},'UniformOutput',false)';

switch usecase
    case 'all'
        % No filtering needed
    case 'all-noshield'
        output(isShield) = [];
    case 'unconv-noshield'
        % Pairs to be used with the unconventional bayes model
        %(i.e. paragliders, gliders,)
        isCostal = contains({listing.folder},{'landuse_beach','landuse_cliff','landuse_volcano','gshhg_gshhs_resf_lvl1'});
        output(~isCostal) = [];
    case 'shield'
        output(~isShield) = [];
    case 'longlinearinfrastructure'
        islli = contains({listing.folder},{'pipeline','roads','waterway','railway','electrictransmission'});
        output(~islli) = [];
    case 'railway'
        irrail = contains({listing.folder},'railway');
        output(~irrail) = [];
    case 'agriculture'
        isAg = contains({listing.folder},{'farm','orchard','vineyard'});
        output(~isAg) = []; 
    otherwise
        error('usecase:unknown','Unknown use case of %s\n',usecase);
end