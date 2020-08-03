% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [sidxFoward, sidxBackward, eidxFoward, eidxBackward, time_s] = findCropIdx(cidx,timeThres_s,method,varargin)

% Todo
% The encounter length CPA check happens outside of this function, it maybe
% worth returning a flag that the time criteria wasn't met

%% Input parser
p = inputParser;

% Required
addRequired(p,'cidx',@(x) isnumeric(x) && x>0); % Index that meets CPA criteria
addRequired(p,'timeThres_s',@(x) isnumeric(x) && x>=0); % How much time to propogate
addRequired(p,'method',@(x) ischar(x) && any(strcmp(x,{'position','time'}))); % Method to use when cropping

% Optional
addOptional(p,'lat',[],@isnumeric); % latitude
addOptional(p,'lon',[],@isnumeric); % longitude
addOptional(p,'v_kts',[],@isnumeric); % velocity
addOptional(p,'time_s',[],@isnumeric); % time

% Parse
parse(p,cidx,timeThres_s,method,varargin{:});

%% Number of input points
switch p.Results.method
    case 'position'
        np = numel(p.Results.lat);
    case 'time'
        np = numel(p.Results.time_s);
end

%% Recalculate time if using velocity or just use input
switch p.Results.method
    case 'position'
        [~, ~,t_legs_s] = calcLegsTime(p.Results.lat,p.Results.lon,p.Results.v_kts);
        time_s = [0; cumsum(t_legs_s)];
    case 'time'
        time_s = p.Results.time_s;
end

%% Determine threshold times relative to time_s(cidx)
sidxTime = time_s(cidx) - timeThres_s;
eidxTime = time_s(cidx) + timeThres_s;

%% Start index
if sidxTime >= 0
    [~,sidxFoward] = min(abs(time_s-sidxTime));
else
    sidxFoward = 1;
end

%% End index
if eidxTime <= max(time_s)
    [~,eidxFoward] = min(abs(time_s-eidxTime));
else
    eidxFoward = np;
end

%% Constant airspeed with fixed spacing which means things are symetric
sidxBackward = eidxFoward;
eidxBackward = sidxFoward;

%% Error checking
assert(sidxFoward <= cidx & cidx <= eidxFoward ,'findCropIdx indexing error')
assert(eidxFoward >= cidx & cidx >= eidxBackward ,'findCropIdx indexing error')

