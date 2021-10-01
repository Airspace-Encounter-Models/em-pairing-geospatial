function [conflicts] = findConflict(mdlType1,mdlType2,track1,track2,lat0_deg,lon0_deg,anchorRange_ft,thresHorz_ft,thresVert_ft,durBeforeConflict_s,durAfterConflict_s,spheroid_ft)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO createEncounters_2

%% Preallocate output
conflicts = double.empty(0,10);

%% Parse
% aircraft 1
ac1Lat_deg = track1.lat_deg;
ac1Lon_deg = track1.lon_deg;
ac1Alt_ft_msl = track1.alt_ft_msl;

% aircraft 2
ac2Lat_deg = track2.lat_deg;
ac2Lon_deg = track2.lon_deg;
ac2Alt_ft_msl = track2.alt_ft_msl;

%% Determine which waypoints are close to anchor point
% Circle a circle around the anchor point
% Find which coordinates are within circle around anchor point
[lat0c,lon0c] = scircle1(lat0_deg,lon0_deg,anchorRange_ft,[],spheroid_ft,'degrees',100);

idxAnchor1 = find(InPolygon(ac1Lat_deg,ac1Lon_deg,lat0c,lon0c));
idxAnchor2 = find(InPolygon(ac2Lat_deg,ac2Lon_deg,lat0c,lon0c));

% When using placeTrack() for positioning Bayes tracks, it is possible to
% have the track originally fly through the anchor point but during
% filtering due to terrain or obstacles, have it removed. This would result
% in the InPolygon() check above returning empty
switch mdlType1
    case 'bayes'
        if isempty(idxAnchor1)
            idxAnchor1 = (1:1:numel(ac1Lat_deg))';
        end
end
switch mdlType2
    case 'bayes'
        if isempty(idxAnchor2)
            idxAnchor2 = (1:1:numel(ac2Lat_deg))';
        end
end

%% Check HMD, VMD, and time criteria
% Make sure both aircraft satisfy the InPolygon check
if ~isempty(idxAnchor1) && ~isempty(idxAnchor2)
    
    % Calculate distance between waypoint combinations
    % Determine which aircraft has more potential points
    if numel(idxAnchor1) >= numel(idxAnchor2)
        % More AC1, so we generate circles for AC2
        idxcol = idxAnchor2;
        latcol = ac2Lat_deg(idxAnchor2);
        loncol = ac2Lon_deg(idxAnchor2);
        
        idxrow = idxAnchor1;
        latrow = ac1Lat_deg(idxAnchor1);
        lonrow = ac1Lon_deg(idxAnchor1);
    else
        % More AC2, so we generate circles for AC1
        idxrow = idxAnchor2;
        latrow = ac2Lat_deg(idxAnchor2);
        lonrow = ac2Lon_deg(idxAnchor2);
        
        idxcol = idxAnchor1;
        latcol = ac1Lat_deg(idxAnchor1);
        loncol = ac1Lon_deg(idxAnchor1);
    end
    
    % Generate scircles for columns
    [latc,lonc] = scircle1(latcol,loncol,repmat(thresHorz_ft,size(loncol)),[],spheroid_ft,'degrees',100);
    
    % Preallocate logical matrix
    isClose = sparse(numel(latrow),numel(latcol));
    
    % Iterate over circles
    for i=1:1:numel(latcol)
        isClose(:,i) = InPolygon(latrow,lonrow,latc(:,i),lonc(:,i));
    end
    
    % Convert from linear indices to subscripts
    [row,col,~] = find(isClose);
    %[row,col] = ind2sub(size(isClose),find(isClose));
    
    % Regenerate c with those that satisfy thresHorz_ft
    % AC1 should be column 1, AC2 should be column 2
    if numel(idxAnchor1) >= numel(idxAnchor2)
        c = [idxrow(row),idxcol(col)];
    else
        c = [idxcol(col),idxrow(row)];
    end
    
    % Calculate HMD, VMD, and filter
    if ~isempty(c)
        % Calculate distance aka HMD
        %hmd_ft = distance(ac1Lat_deg(c(:,1)),ac1Lon_deg(c(:,1)),ac2Lat_deg(c(:,2)),ac2Lon_deg(c(:,2)),spheroid_ft);
        hmd_ft = hypot(track1.north_ft(c(:,1)) - track2.north_ft(c(:,2)),track1.east_ft(c(:,1)) - track2.east_ft(c(:,2)));
        
        % Calculate vertical miss distance
        vmd_ft = abs(ac1Alt_ft_msl(c(:,1)) - ac2Alt_ft_msl(c(:,2)));
        
        % Identify those in VMD and HMD conflict
        lhmd = hmd_ft <= thresHorz_ft;
        lvmd = vmd_ft <= thresVert_ft;
        idxConflicts = find(lhmd & lvmd);
        
        % Filter to those in VMD and HMD conflict
        c = c(idxConflicts,:);
        hmd_ft = hmd_ft(idxConflicts);
        vmd_ft = vmd_ft(idxConflicts);
    end
    
    % Filter time criteria and assign function output
    if ~isempty(c)
        % Total duration of tracks
        len1 = size(track1,1);
        len2 = size(track2,1);
        
        % Identify conflicts that have enough time before
        % and after threshold satisfied
        % This is forward propagation and always valid
        isValid1 = ((c(:,1) + durAfterConflict_s) <= len1) & ((c(:,1) - durBeforeConflict_s) >= 1);
        isValid2 = ((c(:,2) + durAfterConflict_s) <= len2) & ((c(:,2) - durBeforeConflict_s) >= 1);
        
        % If geospatial, also need to account fo backward propagation
        % This could disqualify some valid forward / forward propagation cases
        % additional flexibility can be addressed in a future release
        switch mdlType1
            case 'geospatial'
                isValid1 = isValid1 & ((c(:,1) + durBeforeConflict_s) <= len1) & ((c(:,1) - durAfterConflict_s) >= 1);
        end
        switch mdlType2
            case 'geospatial'
                isValid2 = isValid2 & ((c(:,2) + durBeforeConflict_s) <= len2) & ((c(:,2) - durAfterConflict_s) >= 1);
        end
        
        % Filter
        isValidAll = isValid1 & isValid2;
        c = c(isValidAll,:);
        vmd_ft = vmd_ft(isValidAll);
        hmd_ft = hmd_ft(isValidAll);
        
        % Assign output
        conflicts = zeros(size(c,1),1);
        conflicts(:,1) = ac1Lat_deg(c(:,1)); % AC1 Latitude
        conflicts(:,2) = ac1Lon_deg(c(:,1)); % AC1 Longitude
        conflicts(:,3) = ac2Lat_deg(c(:,2)); % AC2 Latitude
        conflicts(:,4) = ac2Lon_deg(c(:,2)); % AC2 Longitude
        conflicts(:,5) = ac1Alt_ft_msl(c(:,1)); % AC1 Altitude
        conflicts(:,6) = ac2Alt_ft_msl(c(:,2)); % AC2 Altitude
        conflicts(:,7) = c(:,1);
        conflicts(:,8) = c(:,2);
        conflicts(:,9) = hmd_ft;
        conflicts(:,10) = vmd_ft;
    end
end

