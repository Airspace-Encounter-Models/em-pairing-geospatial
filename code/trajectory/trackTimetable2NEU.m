function [arrayOut, speed_ft_s] = trackTimetable2NEU(trackIn,sidx,eidx,dt)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% Converts timetables to arrays in a specific column order
%
% SEE ALSO createEncounters_2 findCropIdx neu2Waypoints timetable

% Input Handling
if nargin < 2; sidx = 1; end;
if nargin < 3; eidx = size(trackIn,1); end
if nargin < 4; dt = seconds(1); end

% Find colummn indicies for NEU variables
switch nargout
    case 1
        idxVar = find(contains(trackIn.Properties.VariableNames,{'north_ft','east_ft','alt_ft_msl'}));
    case 2
        idxVar = find(contains(trackIn.Properties.VariableNames,{'north_ft','east_ft','alt_ft_msl','speed_ft_s'}));
    otherwise
        idxVar = find(contains(trackIn.Properties.VariableNames,{'north_ft','east_ft','alt_ft_msl','speed_ft_s'}));
end

% Filter track timetables
if sidx <= eidx
    ttNEU = trackIn(sidx:1:eidx,idxVar);
else
    ttNEU = trackIn(sidx:-1:eidx,idxVar);
end

% Retime to desired timestep
if ttNEU.Properties.TimeStep ~= dt
    ttNEU = retime(ttNEU,'regular','pchip','TimeStep',dt);
end

% Have filtered NEU track to start at t=0 with positive timestep
ttNEU.Properties.StartTime = seconds(0);

% Convert timetables to arrays in a specific column order
arrayOut = [seconds(ttNEU.Time), ttNEU.north_ft, ttNEU.east_ft, ttNEU.alt_ft_msl];

% Output speed as optional output
if nargout > 1
    speed_ft_s = ttNEU.speed_ft_s;
end
