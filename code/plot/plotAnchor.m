function fig = plotAnchor(anchor,Tdof)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO createEncounters_2

fig = figure;
%geoplot(anchor.lat_deg,anchor.lon_deg,'k*','DisplayName',sprintf('(%0.3f,%0.3f)',anchor.lat_deg,anchor.lon_deg)); hold on;

[uObs,~,ib] = unique(Tdof.obs_type);
for ii=1:1:numel(uObs)
    l = ib == ii;
    [lat,lon] = polyjoin(Tdof.lat_acc_deg(l),Tdof.lon_acc_deg(l)); 
    geoplot(lat,lon,'DisplayName',uObs{ii});hold on;
end
set(gca,'Basemap','satellite','TickLabelFormat','dd');
%legend