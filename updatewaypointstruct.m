% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [waypoints_struct,encCount,encAdd] = updatewaypointstruct(p,waypoints_struct,encCount,ac1NEU_foward,ac2NEU_foward,ac1NEU_backward,ac2NEU_backward)

encAdd = false;

% Adjust timesteps to start at zero
if ~isempty(ac1NEU_foward)
    ac1NEU_foward(:,1) = abs(ac1NEU_foward(:,1) - ac1NEU_foward(1,1));
end
if ~isempty(ac2NEU_foward)
    ac2NEU_foward(:,1) = abs(ac2NEU_foward(:,1) - ac2NEU_foward(1,1));
end
if ~isempty(ac1NEU_backward)
    ac1NEU_backward(:,1) = abs(ac1NEU_backward(:,1) - ac1NEU_backward(1,1));
end
if ~isempty(ac2NEU_backward)
    ac2NEU_backward(:,1) = abs(ac2NEU_backward(:,1) - ac2NEU_backward(1,1));
end

% CSIM requires an update rate of >= 1 second
% Luis: The spacing has to be at least 1 second or greater.
% I know from experience that if you sample at larger time steps
% the trajectories start looking very different,
% but I have found time spacing between 1-3 seconds
% usually does a pretty good job of capturing what you want
[outTime, outX, outY, outZ] = interpTime(p.Results.timeStep_s,ac1NEU_foward(:,1),ac1NEU_foward(:,2),ac1NEU_foward(:,3),ac1NEU_foward(:,4),'isRoundZ',true); ac1NEU_foward = [outTime,outX,outY,outZ];
if ~strcmp(p.Results.acmodel1,'bayes')
    [outTime, outX, outY, outZ] = interpTime(p.Results.timeStep_s,ac1NEU_backward(:,1),ac1NEU_backward(:,2),ac1NEU_backward(:,3),ac1NEU_backward(:,4),'isRoundZ',true); ac1NEU_backward = [outTime,outX,outY,outZ];
end

[outTime, outX, outY, outZ] = interpTime(p.Results.timeStep_s,ac2NEU_foward(:,1),ac2NEU_foward(:,2),ac2NEU_foward(:,3),ac2NEU_foward(:,4),'isRoundZ',true); ac2NEU_foward = [outTime,outX,outY,outZ];
if ~strcmp(p.Results.acmodel2,'bayes')
    [outTime, outX, outY, outZ] = interpTime(p.Results.timeStep_s,ac2NEU_backward(:,1),ac2NEU_backward(:,2),ac2NEU_backward(:,3),ac2NEU_backward(:,4),'isRoundZ',true); ac2NEU_backward = [outTime,outX,outY,outZ];
end

%   WAYPOINTS STRUCTURE:
%   The waypoints structure is an m x n structure matrix, where m is the
%   number of aircraft and n is the number of encounters. This structure
%   matrix has two fields: initial and update. The initial field is a 3
%   element array specifying the north, east, and altitude. The update
%   field is a 4 x n matrix, where n is the number of updates. The rows
%   correspond to the time, north, east, and altitude position of the
%   waypoints.

%%%%%%%%%%%%%%%%%%% COMBO 1
% AC1 forward / AC2 forward
waypoints_ff = neu2wpstruct(ac1NEU_foward,ac2NEU_foward,p.Results.encTime_s,p.Results.initHorz_ft,p.Results.initVert_ft);
if ~isempty(waypoints_ff)
    waypoints_struct = [waypoints_struct waypoints_ff];
    encCount = encCount + 1; encAdd = true;
end

%%%%%%%%%%%%%%%%%%% COMBO 2
% AC1 backward / AC2 foward
waypoints_bf = neu2wpstruct(ac1NEU_backward,ac2NEU_foward,p.Results.encTime_s,p.Results.initHorz_ft,p.Results.initVert_ft);
if ~isempty(waypoints_bf)
    waypoints_struct = [waypoints_struct waypoints_bf];
    encCount = encCount + 1; encAdd = true;
end

% Bayes trajectories can only go foward due to sampling
if ~strcmp(p.Results.acmodel2,'bayes')
    %%%%%%%%%%%%%%%%%%% COMBO 3
    % AC1 foward / AC2 backwards
    waypoints_fb = neu2wpstruct(ac1NEU_foward,ac2NEU_backward,p.Results.encTime_s,p.Results.initHorz_ft,p.Results.initVert_ft);
    if ~isempty(waypoints_fb)
        waypoints_struct = [waypoints_struct waypoints_fb];
        encCount = encCount + 1; encAdd = true;
    end
    
    %%%%%%%%%%%%%%%%%%% COMBO 4
    % AC1 backwards / AC2 backwards
    waypoints_bb = neu2wpstruct(ac1NEU_backward,ac2NEU_backward,p.Results.encTime_s,p.Results.initHorz_ft,p.Results.initVert_ft);
    if ~isempty(waypoints_bb)
        waypoints_struct = [waypoints_struct waypoints_bb];
        encCount = encCount + 1; encAdd = true;
    end
end