function [conflictsOut,sidx1,sidx2,eidx1,eidx2] = findCropIdx(conflicts,mdlType1,mdlType2,track1,track2,initHorz_ft,initVert_ft,durBeforeConflict_s,durAfterConflict_s)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% There are four direction combinations we can arrive at the specific
% waypoint combination. Think of it was foward / backward propagation for
% each aircraft, so 2 X 2 = 4 propagation combinations. This results in
% comprehensive ways to get arrive at this CPA at a given time
%
% SEE ALSO createEncounters_2

%% Input parser
p = inputParser;

% Required
addRequired(p,'conflicts',@isnumeric); %
addRequired(p,'mdlType1',@ischar); %
addRequired(p,'mdlType2',@ischar); %
addRequired(p,'track1',@istimetable); %
addRequired(p,'track2',@istimetable); %
addRequired(p,'initHorz_ft',@isvector); %
addRequired(p,'initVert_ft',@isvector); %
addRequired(p,'durBeforeConflict_s',@isnumeric); %
addRequired(p,'durAfterConflict_s',@isnumeric); %

% Parse
parse(p,conflicts,mdlType1,mdlType2,track1,track2,initHorz_ft,initVert_ft,durBeforeConflict_s,durAfterConflict_s);

%% Preallocate output
sidx1 = double.empty(0,1);
sidx2 = double.empty(0,1);
eidx1 = double.empty(0,1);
eidx2 = double.empty(0,1);

conflictsOut = double.empty(0,10);

%% Forward / Forward Propagation
% Start / end indicies - Forward Propagation - AC1
sidxFwd1 = conflicts(:,7) - durBeforeConflict_s;
eidxFwd1 = conflicts(:,7) + durAfterConflict_s;

% Start / end indicies - Forward Propagation - AC2
sidxFwd2 = conflicts(:,8) - durBeforeConflict_s;
eidxFwd2 = conflicts(:,8) + durAfterConflict_s;

% Potential initial separation - Forward Propagation
r_ft = hypot(track1.east_ft(sidxFwd1)-track2.east_ft(sidxFwd2),track1.north_ft(sidxFwd1)-track2.north_ft(sidxFwd2));
z_ft = abs(track1.alt_ft_msl(sidxFwd1) - track2.alt_ft_msl(sidxFwd2));

% Filter to conflicts that satsify the initial
% separation  conditions with forward
% propagation, as forward propogation will
% always be valid regardless of track type
idx = (r_ft >= initHorz_ft(1)) & (r_ft <= initHorz_ft(2)) & (z_ft >= initVert_ft(1)) & (z_ft <= initVert_ft(2));
idxKeep = find(idx == true);

% Filter
conflictsOut = conflicts(idx,:);
sidxFwd1 = sidxFwd1(idx);
sidxFwd2 = sidxFwd2(idx);
eidxFwd1 = eidxFwd1(idx);
eidxFwd2 = eidxFwd2(idx);

% Assign output
sidx1 = sidxFwd1;
sidx2 = sidxFwd2;
eidx1 = eidxFwd1;
eidx2 = eidxFwd2;

%% Now calculate potential backward propagation
% Start / end indicies - Backward Propagation - AC1
switch mdlType1
    case 'geospatial'
        sidxBck1 = conflictsOut(:,7) + durBeforeConflict_s;
        eidxBck1 = conflictsOut(:,7) - durAfterConflict_s;
    case 'bayes'
        sidxBck1 = double.empty(1,0);
        eidxBck1 = double.empty(1,0);
end

% Start / end indicies - Backward Propagation - AC2
switch mdlType2
    case 'geospatial'
        sidxBck2 = conflictsOut(:,8) + durBeforeConflict_s;
        eidxBck2 = conflictsOut(:,8) - durAfterConflict_s;
    case 'bayes'
        sidxBck2 = double.empty(1,0);
        eidxBck2 = double.empty(1,0);
end

%% Iterate over other propagation combinations
for ii=2:1:4
    % Create indicies pairs
    switch ii
        case 2 % Forward / Backward
            s1 = sidxFwd1; e1 = eidxFwd1;
            s2 = sidxBck2; e2 = eidxBck2;
        case 3 % Backward / Forward
            s1 = sidxBck1; e1 = eidxBck1;
            s2 = sidxFwd2; e2 = eidxFwd2;
        case 4 % Backward / Backward
            s1 = sidxBck1; e1 = eidxBck1;
            s2 = sidxBck2; e2 = eidxBck2;
    end
    
    % Do something if indices are permited
    if ~isempty(s1) && ~isempty(s2)
        % Potential initial separation  - Forward Propagation
        r_ft = hypot(track1.east_ft(s1)-track2.east_ft(s2),track1.north_ft(s1)-track2.north_ft(s2));
        z_ft = abs(track1.alt_ft_msl(s1) - track2.alt_ft_msl(s2));
        
        % Find indices that satify criteria
        idx = (r_ft >= initHorz_ft(1)) & (r_ft <= initHorz_ft(2)) & (z_ft >= initVert_ft(1)) & (z_ft <= initVert_ft(2));
        
        % Assign / append output
        conflictsOut = [conflictsOut;conflictsOut(idx,:)];   
        sidx1 = [sidx1; s1(idx)];
        sidx2 = [sidx2; s2(idx)];
        eidx1 = [eidx1; e1(idx)];
        eidx2 = [eidx2; e2(idx)];
    end
end

