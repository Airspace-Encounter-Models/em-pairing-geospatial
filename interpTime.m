function [outTime, outX, outY, outZ] = interpTime(step_s,time_s,x,y,z,varargin)
% INTERPTIME interpolates X,Y,Z position coordinates with a new timestep
% SEE ALSO: interp1

%% Input parser
p = inputParser;

% Required
addRequired(p,'step_s',@(x) isnumeric(x) && (x >=1)); %  Interpolated timestep, CSIM requires at least >= 1
addRequired(p,'time_s',@isnumeric); %  time
addRequired(p,'x',@isnumeric); %  x position
addRequired(p,'y',@isnumeric); %  y position
addRequired(p,'z',@isnumeric); %  z_position

% Optional
addOptional(p,'isRoundX',false,@islogical); % If true, round outX
addOptional(p,'isRoundY',false,@islogical); % If true, round outY
addOptional(p,'isRoundZ',false,@islogical); % If true, round outZ

% Parse
parse(p,step_s,time_s,x,y,z,varargin{:});

%% Check for unique rows
[~,ia,~] = unique(time_s);

%% Interpolate
% Generate new timesteps
tq_s = (min(time_s(ia)):step_s:max(time_s(ia)))';
outTime = tq_s;

% Find unique rows, because interp1 requires unique vectors
% Due to helicopters, it is potential that position doesn't update
%[~,ia,~] = unique([x,y,z],'rows','stable');

% Interpolate positions
outX = interp1(time_s(ia),x(ia),tq_s);
outY = interp1(time_s(ia),y(ia),tq_s);
outZ = interp1(time_s(ia),z(ia),tq_s);

%% Round
if p.Results.isRoundX; outX = round(outX); end
if p.Results.isRoundY; outY = round(outY); end
if p.Results.isRoundZ; outZ = round(outZ); end

