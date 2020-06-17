function [ac1NEU_foward,ac2NEU_foward,ac1NEU_backward,ac2NEU_backward] = croptracks(p,wypts,encTime_s,ac1Time_s,ac2Time_s,ac1Lat_deg,ac2Lat_deg,ac1Lon_deg,ac2Lon_deg,ac1North_ft,ac2North_ft,ac1East_ft,ac2East_ft,ac1Alt_ft_msl,ac2Alt_ft_msl,altAdjust1_ft,altAdjust2_ft,v1_kts,v2_kts)

% Find CPA  waypoint indices
%cidx1 = find(wypts(1) == ac1Lat_deg & wypts(2) == ac1Lon_deg);
%cidx2 = find(wypts(3) == ac2Lat_deg & wypts(4) == ac2Lon_deg,1,'first');
cidx1 = wypts(7);
cidx2 = wypts(8);

% There are four direction combinations we can arrive at the specific
% waypoint combination. Think of it was foward / backward propagation for
% each aircraft, so 2 X 2 = 4 propagation combinations. This results in
% comprehensive ways to get arrive at this CPA at a given time
switch p.Results.acmodel1
    case 'bayes'
        [sidxFoward1, sidxBackward1,eidxFoward1, eidxBackward1,timesteps1] = findCropIdx(cidx1,encTime_s/2,'time',...
            'time_s',ac1Time_s);
    case 'geospatial'
        [sidxFoward1, sidxBackward1,eidxFoward1, eidxBackward1,timesteps1] = findCropIdx(cidx1,encTime_s/2,'position',...
            'lat',ac1Lat_deg,'lon',ac1Lon_deg,'v_kts',v1_kts);
end

switch p.Results.acmodel2
    case 'bayes'
        [sidxFoward2, sidxBackward2,eidxFoward2, eidxBackward2,timesteps2] = findCropIdx(cidx2,encTime_s/2,'time',...
            'time_s',ac2Time_s);
    case 'geospatial'
        [sidxFoward2, sidxBackward2,eidxFoward2, eidxBackward2,timesteps2] = findCropIdx(cidx2,encTime_s/2,'position',...
            'lat',ac2Lat_deg,'lon',ac2Lon_deg,'v_kts',v2_kts);
end

% AC1
ac1NEU_foward = [timesteps1(sidxFoward1:1:eidxFoward1),ac1North_ft(sidxFoward1:1:eidxFoward1),ac1East_ft(sidxFoward1:1:eidxFoward1),ac1Alt_ft_msl(sidxFoward1:1:eidxFoward1)+altAdjust1_ft];
ac1NEU_backward = [timesteps1(sidxBackward1:-1:eidxBackward1),ac1North_ft(sidxBackward1:-1:eidxBackward1),ac1East_ft(sidxBackward1:-1:eidxBackward1),ac1Alt_ft_msl(sidxBackward1:-1:eidxBackward1)+altAdjust1_ft];

% AC2
ac2NEU_foward = [timesteps2(sidxFoward2:1:eidxFoward2),ac2North_ft(sidxFoward2:1:eidxFoward2),ac2East_ft(sidxFoward2:1:eidxFoward2),ac2Alt_ft_msl(sidxFoward2:1:eidxFoward2)+altAdjust2_ft];
ac2NEU_backward = [timesteps2(sidxBackward2:-1:eidxBackward2),ac2North_ft(sidxBackward2:-1:eidxBackward2),ac2East_ft(sidxBackward2:-1:eidxBackward2),ac2Alt_ft_msl(sidxBackward2:-1:eidxBackward2)+altAdjust2_ft];