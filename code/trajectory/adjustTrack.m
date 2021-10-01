function trackOut = adjustTrack(trackIn, altAdjust_ft, vAdjust_ft_s)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO timetable createEncounters_2 sampleSpeedAlt

%% Input handling
assert(isregular(trackIn),'Input timetable expected to be regular with respect to time, it was not');
assert(trackIn.Properties.TimeStep == seconds(1),'Input timetable timestep expected to be 1 sec, it was not');

%% Preallocate output
trackOut = trackIn;

% Potential variables related to altitude
% Can be multiple because it easier to adjust with the same units (ft)
varsAlt = ["alt_ft_msl", "alt_ft_agl", "up_ft", "down_ft"];

% Variable related to speed
% Only allow one
varSpeed = "speed_ft_s";

% Parse timestep
timeStep_s = seconds(trackIn.Properties.TimeStep);

%% Adjust altitude
if altAdjust_ft ~=0
    % Iterate over potential altitude variables
    for ii = varsAlt
        % Determine if altitude variable ii exists in trackIn
        isField = strcmpi(trackIn.Properties.VariableNames,ii);
        
        % Adjust if variable exists
        if any(isField)
            trackOut.(ii) = trackIn.(ii) + altAdjust_ft;
        end
    end
end

%% Adjust velocity
if vAdjust_ft_s ~=0
    [v_ft_s, d_ft] = estimateSpeedFromGeodetic(trackIn);
    
    % Adjust speed
    trackOut.(varSpeed) = v_ft_s + vAdjust_ft_s;
    
    % Calculate new time required for each leg
    t_legs_s = d_ft ./ abs(trackOut.(varSpeed));
    
    % Account for when v_kts = 0
    t_legs_s(isinf(t_legs_s)) = 0;
    
    % Assume if maneuvering solely vertical, that d_ft = 0,
    % resulting in t_legs_s = NaN
    isVert = isnan(t_legs_s);
    
    % Do something if assumed vertical only maneuver happen
    if any(isVert)
        % Because the timetable is regular and we didn't change vertical rate,
        % we can use the existing timestep
        t_legs_s(isVert) = timeStep_s;
    end
    
    % Update
    trackOut.Time = seconds(cumsum(t_legs_s));
    
    % Remove any duplicate times
    % https://www.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html
    trackOut = sortrows(trackOut);
    trackOut = unique(trackOut);
    dupTimes = sort(trackOut.Time);
    isDup = (diff(dupTimes) == 0);
    if any(isDup)
        trackOut(isDup,:) = [];
    end
    
    % Retime track to have the same regular one second timesteps
    % Also check for potential interpolation errors with retime
    if ~isregular(trackOut)
        trackOut = retime(trackOut,'secondly','pchip');
        isBad = trackOut.lat_deg > 90 | trackOut.alt_ft_msl > 60000 | trackOut.alt_ft_agl > 18000;
        if any(isBad)
            idxLastGood = find(~isBad,1,'last');
            if find(isBad) == size(trackOut,1)
                trackOut{end,:}= trackIn{end,:};
                trackOut.speed_ft_s(end) = trackOut.speed_ft_s(idxLastGood);
            else
                trackOut = trackOut(1:idxLastGood,:);
            end
        end
    end
    
        % Have track start at t = 0
    if trackOut.Properties.StartTime ~= seconds(0)
        trackOut.Properties.StartTime = seconds(0);
    end
    
end