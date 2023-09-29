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
idx_iso = find(contains({ne_admin.iso_3166_2},iso_3166_2,'IgnoreCase',true));
if isempty(idx_iso)
    error('iso_3166_2:missing','Cannot find iso_3166_2 = %s in natural earth data', iso_3166_2);
end
BoundingBox_wgs84 = ne_admin(idx_iso).BoundingBox;

neLat_deg = ne_admin(idx_iso).Lat;
neLon_deg = ne_admin(idx_iso).Lon;

% Create grid
[gridLat_deg, gridLon_deg] = meshgrid(min(neLat_deg):nm2deg(anchorRange_nm):max(neLat_deg),min(neLon_deg):nm2deg(anchorRange_nm):max(neLon_deg));

% Filter grid for those only in the location
isGrid = inpolygon(gridLat_deg,gridLon_deg,neLat_deg,neLon_deg);

gridLat_deg = gridLat_deg(isGrid);
gridLon_deg = gridLon_deg(isGrid);

if isPlot
   figure; set(gcf,'name',['createISO31662Grid: ' iso_3166_2]);
   
   geoscatter(gridLat_deg,gridLon_deg,[],categorical(cellstr(class)),'.');
   set(gca,'Basemap','streets-dark','TickLabelFormat','dd');
end