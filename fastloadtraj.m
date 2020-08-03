% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [time_s, lat_deg, lon_deg, alt_ft_msl, alt_ft_agl] = textscantraj(inFile)

% Open, read, and close file
% We ignore the last column with the airspace name because of a bug
% where some airspace names have a comma and MATLAB thinks that is
% another delimiter
fid = fopen(inFile,'r');
textRaw = textscan(fid,'%f %f %f %f %f %f %f %f %f %s %*[^\n]','EndOfLine','\n','Delimiter',',','HeaderLines',1);
fclose(fid);

% Parse
% Header order assumed to be:
% time_s,heading_deg,groundspeed_kt,el_ft_msl,alt_ft_msl,alt_ft_agl,altRate_fps,lat_deg,lon_deg,AirspaceClass,AirspaceName
time_s = textRaw{1};
lat_deg = textRaw{8};
lon_deg = textRaw{9};
alt_ft_msl = textRaw{5};
alt_ft_agl = textRaw{6};