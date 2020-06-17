function [v1_kts,v2_kts,altAdjust1_ft,altAdjust2_ft] = samplespeedalt(p,n,ac1Alt_ft_agl,ac2Alt_ft_agl)

%% Sample aircraft speed for each cpa combination
% This should eventually be for each timestep, not a fixed value
switch p.Results.acmodel1
    case 'bayes'
        v1_kts = nan(n,1);
    case 'geospatial'
        v1_kts = randi(p.Results.vrange1_kts,n,1);
end
switch p.Results.acmodel2
    case 'bayes'
        v2_kts = nan(n,1);
    case 'geospatial'
        v2_kts = randi(p.Results.vrange1_kts,n,1);
end

%% Sample altitude adjustment for each cpa combination (geospatial only)
% Get AGL altitude
% Calculate the interval of potential altitude adjustments
% Sample the interval and set altitude adjustments for each combination
if p.Results.isSampleAlt1
    below = p.Results.minAlt1_ft_agl-min(ac1Alt_ft_agl);
    above = p.Results.maxAlt1_ft_agl-max(ac1Alt_ft_agl);
    if below <= above
        adjustSpan_ft_agl = below:25:above;
    else
        adjustSpan_ft_agl = above:-25:(below + above); % Need to bring aircraft down
    end
    altAdjust1_ft = adjustSpan_ft_agl(randi(numel(adjustSpan_ft_agl),n,1))';
else
    altAdjust1_ft = zeros(n,1);
end

if p.Results.isSampleAlt2
    below = p.Results.minAlt2_ft_agl-min(ac2Alt_ft_agl);
    above = p.Results.maxAlt2_ft_agl-max(ac2Alt_ft_agl);
    if below <= above
        adjustSpan_ft_agl = below:25:(below + above);
    else
        adjustSpan_ft_agl = above:-25:below; % Need to bring aircraft down
    end
    altAdjust2_ft = adjustSpan_ft_agl(randi(numel(adjustSpan_ft_agl),n,1))';
else
    altAdjust2_ft = zeros(n,1);
end

%% Make sure outputs are columns
if ~iscolumn(v1_kts); v1_kts = v1_kts'; end;
if ~iscolumn(v2_kts); v2_kts = v2_kts'; end;
if ~iscolumn(altAdjust1_ft); altAdjust1_ft = altAdjust1_ft'; end;
if ~iscolumn(altAdjust2_ft); altAdjust2_ft = altAdjust2_ft'; end;