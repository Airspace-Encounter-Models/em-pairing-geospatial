function [time_s,lat_deg,lon_deg,alt_ft_msl,alt_ft_agl,xNorth_ft,yEast_ft,zDown_ft,isObstacle] = placeTrack(inFile,lat0,lon0,varargin)


%% Set up input parser
p = inputParser;

% Required
addRequired(p,'inFile');
addRequired(p,'lat0',@(x) isnumeric(x) && numel(x) == 1);
addRequired(p,'lon0',@(x) isnumeric(x) && numel(x) == 1);

% Optional
addOptional(p,'spheroid',wgs84Ellipsoid('ft'));
addOptional(p,'Delimiter',',');
addOptional(p,'z_units','agl',@(x) isstr(x) && any(strcmpi(x,{'agl','msl'})));
addOptional(p,'z_agl_tol_ft',250,@(x) isnumeric(x) && numel(x) == 1);

% Optional - DEM
addOptional(p,'dem','globe',@ischar);
addOptional(p,'Z_m',{},@isnumeric);
addOptional(p,'refvec',{},@isnumeric);

% Optional - Obstacles
addOptional(p,'latObstacle',[],@isnumeric);
addOptional(p,'lonObstacle',[],@isnumeric);
addOptional(p,'altObstacle_ft_agl',[],@isnumeric);

% Optional - inFile parsing
addOptional(p,'labelTime','time_s');
addOptional(p,'labelX','x_ft');
addOptional(p,'labelY','y_ft');
addOptional(p,'labelZ','z_ft');

% Optional Plot
addOptional(p,'isPlot',false,@islogical);

% Parse
parse(p,inFile,lat0,lon0,varargin{:});

%% Load trajectory
T = readtable(inFile,'Delimiter',p.Results.Delimiter);

% Interpolate
[tq_s, x_ft, y_ft, z_ft] = interpTime(1,T.(p.Results.labelTime),T.(p.Results.labelX),T.(p.Results.labelY),T.(p.Results.labelZ));
time_s = tq_s;
zDown_ft = z_ft;

% Preallocate, although this will be overwritten
alt_ft_agl = z_ft;

% Velocity
%v2_kts = interp1(T.time_s,T.groundspeed_kt,tq_s);
%if any(v2_kts == 0); v2_kts(v2_kts == 0) = 1e-3; end

%% Filter and parse out obstacles
l = p.Results.altObstacle_ft_agl >= min(zDown_ft);
latObstacle = p.Results.latObstacle(:,l);
lonObstacle = p.Results.lonObstacle(:,l);
numObstacle = size(latObstacle,2);

%% Translate and rotate trajectory so it doesn't hit an obstacle
c = 1;
isObstacle = true; % Set true to enter while loop
while c <= 10 & isObstacle
    % Select which index will be (0,0)
    % Whatever is (0,0) will pass directly through the anchorpoint
    i0 = randi(numel(time_s),1,1);
    x_ft = x_ft - x_ft(i0);
    y_ft = y_ft - y_ft(i0);
    
    % Randomly rotate and assign
    % https://en.wikipedia.org/wiki/Rotation_matrix
    rDeg = randi(360,1,1);
    yEast_ft = x_ft * cosd(rDeg) - y_ft * sind(rDeg);
    xNorth_ft = x_ft * sind(rDeg) + y_ft * cosd(rDeg);
    
    % Convert to lat / lon
    % h(i0) approx el_ft_msl - z_ft(i0);
    h0 = 0; % placeholder, not really used
    [lat_deg,lon_deg,~] = ned2geodetic(x_ft,y_ft,z_ft,lat0,lon0,h0,p.Results.spheroid);
    
    % Check to make sure not hitting obstacle
    isObstacle = false(numObstacle,1);
    for i=1:1:numObstacle
        isObstacle(i) = any(InPolygon(lat_deg,lon_deg,latObstacle(:,i),lonObstacle(:,i)));
    end
    isObstacle = any(isObstacle);
    
    % There is a horizontal conflict, so we need to check vertical
    if isObstacle
        isObstacle = any(any(p.Results.altObstacle_ft_agl(isObstacle)' >= zDown_ft));
    end
    
    % Advance counter
    c = c + 1;
end

if isObstacle
    warning('placeTracks:obstacle','After %i tries, Bayes track still hitting obstacles\nFILE = %s\n',10,inFile);
end

%% Convert to MSL if needed and assign
alt_ft_msl = nan(size(lat_deg));
switch p.Results.z_units
    case 'agl'
        % Get elevation lat / lon
        if ~any(strcmpi(p.UsingDefaults,'Z_m'))
            [el_ft_msl,~,~,~] = msl2agl(lat_deg, lon_deg,p.Results.dem,'Z_m',p.Results.Z_m,'refvec',p.Results.refvec);
        else
            [el_ft_msl,~,~,~] = msl2agl(lat_deg, lon_deg, p.Results.dem);
        end
        
        % msl2agl() may return an empty array if there are too many NaN in
        % the dem, if this happens, throw everything out
        if ~isempty(el_ft_msl)
            % Want to persve the the correct AGL at the anchor point
            %[~,I] = min(z_ft);
            I = i0;
            
            % Convert from AGL to MSL
            % Bayes uncorrelated & nonconventional model use AGL units in the lowest
            % altitude bins, so this calculation perseves the AGL shape
            % However Bayes HAA uses MSL, so we're forcing a MSL = AGL assummption for HAA
            alt_ft_msl = z_ft + el_ft_msl(I);
            
            % Determine which points are underground and bad
            alt_ft_agl = alt_ft_msl - el_ft_msl;
            isGood = alt_ft_agl >=0 & abs(alt_ft_agl-z_ft) <= p.Results.z_agl_tol_ft;
        else
            isGood = false(size(time_s));
        end
        
        % Find longest sequence of consecutive non-zero values
        % https://www.mathworks.com/matlabcentral/answers/404502-find-the-longest-sequence-of-consecutive-non-zero-values-in-a-vector#answer_323627
        zpos = find(~[0 isGood' 0]);
        [~, grpidx] = max(diff(zpos));
        idx = zpos(grpidx):zpos(grpidx+1)-2;
        
        % Filter using idx
        % This produces the longest track above ground
        time_s = time_s(idx);
        lat_deg = lat_deg(idx);
        lon_deg = lon_deg(idx);
        alt_ft_msl = alt_ft_msl(idx);
        alt_ft_agl = alt_ft_agl(idx);
        xNorth_ft = xNorth_ft(idx);
        yEast_ft = yEast_ft(idx);
        zDown_ft = zDown_ft(idx);
        
    case 'msl'
        % z_ft is already in MSL
        alt_ft_agl = z_ft;
        alt_ft_msl = z_ft;
        isGood = true(size(z_ft));
end

%% Plot
if p.Results.isPlot
    figure; set(gcf,'name',inFile);
    subplot(1,2,1);
    geoshow(lat_deg,lon_deg); grid on;
    xlabel('Latitude (deg)'); ylabel('Longitude (deg)');
    
    subplot(1,2,2);
    plot(tq_s,z_ft,tq_s,alt_ft_agl,tq_s,el_ft_msl,time_s,alt_ft_msl); grid on;
    xlabel('time (s)');
    legend('z_ft','alt_ft_agl','el_ft_msl','alt_ft_msl','Interpreter','none','Location','best');
end

%% Make sure time starts at 1
if any(isGood) && time_s(1) ~=1
    time_s = (time_s - time_s(1))+1;
end
