function [gridLat_deg, gridLon_deg, BoundingBox_wgs84] = createISO31662Grid(iso_3166_2, anchorRange_nm, isPlot)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO preallocAnchors

% Input Handling
if nargin < 3; isPlot = false; end

% Load Natural Earth Adminstrative Boundaries
ne_admin = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'NE-Adminstrative' filesep 'ne_10m_admin_1_states_provinces.shp'],'UseGeoCoords',true);

% Create bounding box
BoundingBox_wgs84 = ne_admin(strcmp({ne_admin.iso_3166_2},iso_3166_2)).BoundingBox;

neLat_deg = ne_admin(strcmp({ne_admin.iso_3166_2},iso_3166_2)).Lat;
neLon_deg = ne_admin(strcmp({ne_admin.iso_3166_2},iso_3166_2)).Lon;

% Create grid
[gridLat_deg, gridLon_deg] = meshgrid(min(neLat_deg):nm2deg(anchorRange_nm):max(neLat_deg),min(neLon_deg):nm2deg(anchorRange_nm):max(neLon_deg));

% Filter grid for those only in the location
isGrid = InPolygon(gridLat_deg,gridLon_deg,neLat_deg,neLon_deg);

gridLat_deg = gridLat_deg(isGrid);
gridLon_deg = gridLon_deg(isGrid);

if isPlot
   figure; set(gcf,'name',['createISO31662Grid: ' iso_3166_2]);
   
   geoscatter(gridLat_deg,gridLon_deg,[],categorical(cellstr(class)),'.');
   set(gca,'Basemap','streets-dark','TickLabelFormat','dd');
end