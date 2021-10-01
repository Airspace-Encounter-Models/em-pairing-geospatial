function track = loadTrack(mdlType, arg2, isStartEndNonZero, rangeAlt_ft_agl, rangeV_ft_s, spheroid_ft, lat0_deg,lon0_deg,el0_ft_msl,Z_m,refvec,R,Tdof,dofMaxRange_ft,dofMaxVert_ft,bayesDuration_s,z_agl_tol_ft)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO timetable textscantraj UncorEncounterModel/track

% Inputs hardcode
% Variable name associated with speed
varSpeed = 'speed_ft_s';

% Behavior depends on model type
switch mdlType
    case 'geospatial'
        % Input handling, expected to be a filename
        assert(isfile(arg2));
        
        % Load track
        track = textscantraj(arg2);
        track.Properties.StartTime = seconds(0);
        
        % em-model-geospatial currently does not calculate east, north, up
        % So calculate east, north, up here
        [track.east_ft,track.north_ft,~] = geodetic2enu(track.lat_deg,track.lon_deg,track.alt_ft_agl,lat0_deg,lon0_deg,el0_ft_msl,spheroid_ft);
        isAltGood = true;
        
    case 'bayes'
        % Input handling, expected to be a class object
        assert(isa(arg2,'UncorEncounterModel'));
        
        isAltGood = false;
        c = 1;
        while ~isAltGood && c < 10
            
            track = arg2.track(1,bayesDuration_s,'lat0_deg',lat0_deg,'lon0_deg',lon0_deg,'coordSys','geodetic',...
                'Tdof',Tdof,'dofMaxRange_ft',dofMaxRange_ft,'dofMaxVert_ft',dofMaxVert_ft,'Z_m',Z_m,'refvec',refvec,'R',R,'z_agl_tol_ft',z_agl_tol_ft);
            track = track{1};
            
            % Check that altitude range is satsified for Bayes
            if ~isempty(track) && (min(track.alt_ft_agl) >= rangeAlt_ft_agl(1) && max(track.alt_ft_agl) <= rangeAlt_ft_agl(2))
                isAltGood = true;
            else
                c = c + 1;
            end
        end
end

if ~isAltGood
    track = timetable.empty(0,12);
else
    
    % Remove points with zero altitude at start and end
    % This often due to a vertical climb / descend at the start or end
    if isStartEndNonZero
        isZero = floor(track.('alt_ft_agl')) == 0;
        sidx = find(~isZero,1,'first');
        eidx = find(~isZero,1,'last');
        track = track(sidx:eidx,:);
    end
    
    % Remove any duplicate times
    % https://www.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html
    track = sortrows(track);
    track = unique(track);
    dupTimes = sort(track.Time);
    isDup = (diff(dupTimes) == 0);
    if any(isDup)
        track(isDup,:) = [];
    end
    
    % Retime tracks to have the same regular one second timesteps
    track = retime(track,'secondly','pchip');
    
    % Calculate speed
    switch mdlType
        case 'geospatial'
            % If not available, estimate it from position
            if ~any(strcmpi(varSpeed,track.Properties.VariableNames))
                track.(varSpeed) = estimateSpeedFromGeodetic(track);
            end
    end
    
    % Remove points with zero speed at start and end
    % This often due to a vertical climb / descend at the start or end
    if isStartEndNonZero
        isZero = floor(track.(varSpeed)) == 0;
        sidx = find(~isZero,1,'first');
        eidx = find(~isZero,1,'last');
        track = track(sidx:eidx,:);
    end
    
    % Have track start at t = 0
    if track.Properties.StartTime ~= seconds(0)
        track.Properties.StartTime = seconds(0);
    end
end