function [altAdjust_ft, vAdjust_ft_s] = sampleSpeedAlt(mdlType,track,rangeAlt_ft_agl,rangeV_ft_s, isSampleAlt, isSampleV,n)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Input Handling
% Never assume we want adjustment, user must be explicit
if nargin < 5; isSampleAlt = false; end
if nargin < 6; isSampleV = false; end

%% Sample aircraft speed for n times
% Behavior depends on model type
switch mdlType
    case 'geospatial'
        if isSampleV
            imin = min([0 rangeV_ft_s(1) - min(track.speed_ft_s)]);
            imax = max([0 rangeV_ft_s(2) - max(track.speed_ft_s)]);
            
            imin = ceil(imin);
            imax = floor(imax);
            
            if imax < imin; imax = imin; end    
            vAdjust_ft_s = randi([imin imax],n,1);
        else
            vAdjust_ft_s = zeros(n,1);
        end
    case 'bayes'
        vAdjust_ft_s = zeros(n,1);
end

%% Sample altitude adjustment n times (geospatial only)

% Behavior depends on model type
switch mdlType
    case 'geospatial'
        % Get AGL altitude
        % Calculate the interval of potential altitude adjustments
        % Sample the interval and set altitude adjustments for each combination
        if isSampleAlt
            below = rangeAlt_ft_agl(1)-min(track.alt_ft_agl);
            above = rangeAlt_ft_agl(2)-max(track.alt_ft_agl);
            if below <= above
                adjustSpan_ft_agl = below:25:above;
            else
                adjustSpan_ft_agl = above:-25:(below + above); % Need to bring aircraft down
            end
            
            if ~isempty(adjustSpan_ft_agl)
                altAdjust_ft = adjustSpan_ft_agl(randi(numel(adjustSpan_ft_agl),n,1))';
            else
                altAdjust_ft = zeros(n,1);
            end
        else
            altAdjust_ft = zeros(n,1);
        end
        
    case 'bayes'
        altAdjust_ft = zeros(n,1);
end

%% Make sure outputs are columns
if ~iscolumn(altAdjust_ft); altAdjust_ft = altAdjust_ft'; end
if ~iscolumn(vAdjust_ft_s); vAdjust_ft_s = vAdjust_ft_s'; end
