function [v_ft_s, d_ft] = estimateSpeedFromGeodetic(track)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

% Parse timestep
timeStep_s = seconds(track.Properties.TimeStep);

% Calculate distance between legs
[~,d_nm] = legs(track.lat_deg,track.lon_deg,'rh');

% Convert from nautical miles to feet
d_ft = d_nm * 6076.12;
d_ft = [d_ft(1); d_ft];
d_ft = abs(d_ft);

% Estimate speed
v_ft_s =  d_ft./ timeStep_s;

