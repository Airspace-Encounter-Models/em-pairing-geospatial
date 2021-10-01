function [anchors, gridLat_deg, gridLon_deg, airspace] = preallocAnchors(iso_3166_2, anchorRange_nm, fileAirspace, isPlot)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO createISO31662Grid IdentifyGeographicVariable

% Input handling
if nargin < 3 || isempty(fileAirspace)
    fileAirspace = [getenv('AEM_DIR_CORE') filesep 'output' filesep 'airspace.mat'];
end
if nargin < 4; isPlot = false; end

% Create grid
[gridLat_deg, gridLon_deg, BoundingBox_wgs84] = createISO31662Grid(iso_3166_2, anchorRange_nm, false);

% Calculate geographic domain
G = IdentifyGeographicVariable(gridLat_deg,gridLon_deg);

% Load airspace
load(fileAirspace,'airspace');

% Filter to include airspace only in bounding box
[~, ~, inAir] = filterboundingbox(airspace.LAT_deg,airspace.LON_deg,BoundingBox_wgs84);

% Filter airspace based on altitude surface
isAir = inAir & cellfun(@min,airspace.LOWALT_ft_agl) <= 0;
airspace = airspace(isAir,:);

% Iterate over airspace
class = repmat('O',size(gridLon_deg));
for  ii=1:1:size(airspace,1)
    isAirspace = InPolygon(gridLat_deg,gridLon_deg,airspace.LAT_deg{ii},airspace.LON_deg{ii});
    class(isAirspace) = airspace.CLASS(ii);
end

% Plot
if isPlot
    figure; set(gcf,'name',['preallocAnchors: ' iso_3166_2]);
    isB = class == 'B';
    isC = class == 'C';
    isD = class == 'D';
    isO = class == 'O';
    
    geoscatter(gridLat_deg(isB),gridLon_deg(isB),'r.','DisplayName','Class B'); hold on;
    geoscatter(gridLat_deg(isC),gridLon_deg(isC),'m.','DisplayName','Class C');
    geoscatter(gridLat_deg(isD),gridLon_deg(isD),'b.','DisplayName','Class D');
    geoscatter(gridLat_deg(isO),gridLon_deg(isO),'w.','DisplayName','Other'); hold off;
    set(gca,'Basemap','streets-dark','TickLabelFormat','dd','FontSize',12,'FontName','Arial','FontWeight','bold');
    legend('Location','best','TextColor','w','FontSize',10,'FontName','Arial','FontWeight','bold'); legend('boxoff');
end

% Assign output
anchors = table(gridLat_deg,gridLon_deg,class,cell(size(gridLon_deg)),zeros(size(gridLon_deg)),G,'VariableNames',{'lat_deg','lon_deg','class','files','num_geospatial','G'});
