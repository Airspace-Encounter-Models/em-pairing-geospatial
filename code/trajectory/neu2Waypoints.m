function wp = neu2Waypoints(ac1,ac2)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%   WAYPOINTS STRUCTURE:
%   The waypoints structure is an m x n structure matrix, where m is the
%   number of aircraft and n is the number of encounters. This structure
%   matrix has two fields: initial and update. The initial field is a 3
%   element array specifying the north, east, and altitude. The update
%   field is a 4 x n matrix, where n is the number of updates. The rows
%   correspond to the time, north, east, and altitude position of the
%   waypoints.

% Modify time so it starts from 0
if ac1(1,1) ~= 0; ac1(:,1) = ac1(:,1) - ac1(1,1); end
if ac2(1,1) ~= 0; ac2(:,1) = ac2(:,1) - ac2(1,1); end

% Fill in wp struct for ac1
for idx = 1:size(ac1,1)
    if idx == 1
        wp(1,1).initial = [ac1(idx,2);...
            ac1(idx,3);...
            ac1(idx,4);...
            ];
    else
        wp(1,1).update(:,idx-1) = [ac1(idx,1);...
            ac1(idx,2);...
            ac1(idx,3);...
            ac1(idx,4);...
            ];
    end
end

% Fill in wp struct for ac2
for idx = 1:size(ac2,1)
    if idx == 1
        wp(2,1).initial = [ac2(idx,2);...
            ac2(idx,3);...
            ac2(idx,4);...
            ];
    else
        wp(2,1).update(:,idx-1) = [ac2(idx,1);...
            ac2(idx,2);...
            ac2(idx,3);...
            ac2(idx,4);...
            ];
    end
end
