% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function waypoints_csim = neu2wpstruct(ac1NEU,ac2NEU,encTime_s,initHorz_ft,initVert_ft)

waypoints_csim = [];

% First check encounter duration
if ac1NEU(end,1) >= encTime_s && ac2NEU(end,1) >= encTime_s
    
    % Calculate
    range_ft = sqrt((ac1NEU(1,2)-ac2NEU(1,2))^2 + (ac1NEU(1,3)-ac2NEU(1,3))^2);
    z_ft = abs(ac1NEU(1,4) - ac2NEU(1,4));
    
    % Evaluate range and vertical sep.
    if initHorz_ft(1) <= range_ft&&...
            initHorz_ft(2) >= range_ft&&...
            initVert_ft(1) <= z_ft &&...
            initVert_ft(2) >= z_ft
        
        waypoints_csim = neu_to_wp_struct(ac1NEU,ac2NEU);
    else
        %fprintf('Skipping because initial separation conditions not met when i=%i, j=%i\n',i,j);
    end
end