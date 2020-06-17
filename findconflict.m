function [wypts,ac1Time_s,ac1Lat_deg,ac1Lon_deg,ac1Alt_ft_msl,ac2Time_s,ac2Lat_deg,ac2Lon_deg,ac2Alt_ft_msl] = findconflict(p,time1_s,lat1_deg,lon1_deg,alt1_ft_msl,time2_s,lat2_deg,lon2_deg,alt2_ft_msl,lat0,lon0,anchorRange_ft)

%% Get waypoints
% aircraft 1
ac1Time_s = time1_s(p.Results.rowstart1:end-p.Results.rowend1);
ac1Lat_deg = lat1_deg(p.Results.rowstart1:end-p.Results.rowend1);
ac1Lon_deg = lon1_deg(p.Results.rowstart1:end-p.Results.rowend1); %#ok
ac1Alt_ft_msl = alt1_ft_msl(p.Results.rowstart1:end-p.Results.rowend1); %#ok

% aircraft 2
ac2Time_s = time2_s(p.Results.rowstart2:end-p.Results.rowend2); %#ok
ac2Lat_deg = lat2_deg(p.Results.rowstart2:end-p.Results.rowend2); %#ok
ac2Lon_deg = lon2_deg(p.Results.rowstart2:end-p.Results.rowend2); %#ok
ac2Alt_ft_msl = alt2_ft_msl(p.Results.rowstart2:end-p.Results.rowend2); %#ok

%% Determine which waypoints are close to anchor point
% Circle a circle around the anchor point
% Find which coordinates are within circle around anchor point
[lat0c,lon0c] = scircle1(lat0,lon0,anchorRange_ft,[],p.Results.spheroid,'degrees',100);

idxAnchor1 = find(InPolygon(ac1Lat_deg,ac1Lon_deg,lat0c,lon0c));
idxAnchor2 = find(InPolygon(ac2Lat_deg,ac2Lon_deg,lat0c,lon0c));

% When using placeTrack() for positioning Bayes tracks, it is possible to
% have the track originally fly through the anchor point but during
% filtering due to terrain or obstacles, have it removed. This would result
% in the InPolygon() check above returning empty
switch p.Results.acmodel1
    case 'bayes'
        if isempty(idxAnchor1)
            idxAnchor1 = (1:1:numel(ac1Lat_deg))';
        end
end
switch p.Results.acmodel2
    case 'bayes'
        if isempty(idxAnchor2)
            idxAnchor2 = (1:1:numel(ac2Lat_deg))';
        end
end

% Deprecated using distance
% This way is (negligible more) accurate because it takes into account the spheriod but is also much slower
%idxAnchor1 = find(distance(ac1Lat_deg,ac1Lon_deg,lat0,lon0,p.Results.spheroid) <= anchorRange_ft);
%idxAnchor2 = find(distance(ac2Lat_deg,ac2Lon_deg,lat0,lon0,p.Results.spheroid) <= anchorRange_ft);

% DEPRECATED - REMOVE IN FUTURE UPDATES
% Calculate combinations of waypoints
% Filter combinations to those only near the anchor point
% c = combvec(idxAnchor1',idxAnchor2')';
%
% % Randomly order combinations
% % Downsample combination of waypoints if needed
% if size(c,1) > p.Results.maxCombo
%     pidx = randperm(size(c,1),p.Results.maxCombo);
%     fprintf('Too many waypoints combinations for pair %i, sampling %i combinations\n',i,p.Results.maxCombo);
% else
%     pidx = randperm(size(c,1),size(c,1));
% end
% c = c(pidx',:);

%% Return if not both aircraft satisfy the InPolygon check
if isempty(idxAnchor1) | isempty(idxAnchor2)
    c = zeros(0,2);
    wypts = zeros(0,8);
    return;
end

%% Calculate distance between waypoint combinations
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
[latc,lonc] = scircle1(latcol,loncol,repmat(p.Results.thresHorz_ft,size(loncol)),[],p.Results.spheroid,'degrees',100);

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

% Filter to desired number of encounters
if size(c,1) > p.Results.maxEncPerPair
    c = c(randperm(size(c,1),p.Results.maxEncPerPair),:);
end

% Calculate distance (uncomment for debugging)
% d_ft = distance(ac1Lat_deg(c(:,1)),ac1Lon_deg(c(:,1)),ac2Lat_deg(c(:,2)),ac2Lon_deg(c(:,2)),p.Results.spheroid);

%% Organize waypoints
wypts = zeros(size(c,1),8);
if ~isempty(wypts)
    wypts(:,1) = ac1Lat_deg(c(:,1)); % AC1 Latitude
    wypts(:,2) = ac1Lon_deg(c(:,1)); % AC1 Longitude
    wypts(:,3) = ac2Lat_deg(c(:,2)); % AC2 Latitude
    wypts(:,4) = ac2Lon_deg(c(:,2)); % AC2 Longitude
    wypts(:,5) = ac1Alt_ft_msl(c(:,1)); % AC1 Altitude
    wypts(:,6) = ac2Alt_ft_msl(c(:,2)); % AC2 Altitude
    wypts(:,7) = c(:,1);
    wypts(:,8) = c(:,2);
    % wypts(:,9) = d_ft;
end