% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function suas_traj_to_wp(inputDir, fout, varargin)
    if ~contains(fout, '.dat')
        error('Output should be .dat file');
    end
    opts = inputParser;
    opts.addParameter('numTrialsStart', 3e4, @isnumeric);
    opts.addParameter('numTrialsEnd', 1e4, @isnumeric);
    opts.addParameter('minSpeedOwn',25, @isnumeric); % kts; CDF resampled until speed is >= minSpeed
    opts.addParameter('minSpeedInt',0, @isnumeric);
    opts.addParameter('tcpa',80, @isnumeric); % time of cpa 
    opts.addParameter('speedDistributionFile','suasJcruiseCDF.mat',@ischar); %.mat file containing sUAS speed distributions
    opts.addParameter('runPar',true, @isboolean);
    opts.addParameter('minAlt',1200, @isnumeric);
    opts.addParameter('speed',60, @isnumeric); %in knots
    opts.addParameter('minT', 80, @isnumeric); % minimum time of encounters desired
    opts.addParameter('maxT',120, @isnumeric); % maximum time of encounters desired
    opts.addParameter('tol',5,@isnumeric) % tolerance in seconds for removing encounters with tca not within tol
    opts.parse(varargin{:});
    
    numTrialsStart = opts.Results.numTrialsStart;
    numTrialsEnd = opts.Results.numTrialsEnd;
    minSpeedOwn = opts.Results.minSpeedOwn;
    minSpeedInt = opts.Results.minSpeedInt;
    tcpa = opts.Results.tcpa;
    speedDistributionFile = opts.Results.speedDistributionFile;
    if ~exist(speedDistributionFile, 'file')
        speedDistributionFile = fullfile(mitcas_paths.mitcas_shared, 'encounter_files','sUAS','source','data_files','suasJcruiseCDF.mat');
    end
    runPar = opts.Results.runPar;
    minAlt = opts.Results.minAlt;
    speed = opts.Results.speed;
    minT = opts.Results.minT;
    maxT = opts.Results.maxT;
    tol = opts.Results.tol;
    
    temp = strsplit(fout,'.');
    fout = temp{1};
    
    %% create encounters
    fprintf('Creating encounters... \n');
    tic
    createEncounters(inputDir, fout, speed, minT, maxT);
    toc
    
    fout = [fout '_WP'];
    %% pairs suas trajectories with each other
    fprintf('Pairing sUAS trajectories...\n');
    tic
    pair_uas2_wp([fout '.dat'], [fout 'Paired.dat'],numTrialsStart, minSpeedOwn, minSpeedInt, tcpa, speedDistributionFile, runPar);
    toc
    
    %% adjust waypoints
    fprintf('Adjusting waypoints... \n');
    tic
    adjWP([fout 'Paired.dat'], [fout 'Adj.dat'], tcpa, minAlt);
    toc
    
    %% finalize encounter sets
    fprintf('Finalizing encounter sets... \n');
    tic
    finWP([fout 'Adj.dat'], [fout 'Final.dat'], numTrialsEnd, tcpa, tol);
    toc
    
    fprintf('Encounter set finished!\n');
    
end

function createEncounters(dirName, fout, speed_kts, minT, maxT)
    % This script takes data from .csv files and creates trajectories

    % Inputs
    speed = speed_kts * 1.68781;
    R = 3959*5280; % radius of earth, in feet (assuming perfect sphere)
    foutWP = [fout '_WP.dat'];
    foutS = [fout '_S.dat'];

    d = dir(dirName);
    fNum = sum(~[d.isdir]); % number of files in directory
    filenames = {};
    for ii = 1:numel(d)
        if contains(d(ii).name,'.csv')
            filenames = [filenames fullfile(d(ii).folder,d(ii).name)];
        end
    end

    numAc = 2;
    numEnc = fNum; % set to fNum to use all files
    if numel(filenames) ~= numEnc
        error('Number of filenames doesn''t match numEnc in createEncounters method');
    end
    waypoints = struct('initial', cell(numAc, 1), 'update', cell(numAc, 1));
    index = 1; % array index increments separately from main loop

    for i = 1:numel(filenames)
        fname = filenames{i};
%         if exist(fname,'file')~=2 % file doesn't exist
%             fprintf('File number %g does not exist.\n', i);
%             continue;
%         end
        % Read raw input data
        data = csvread(fname, 1, 0); % skips header row
        % Remove climb and descend stages
        data = data(2:end-1, :);
        % Remove duplicate time steps
        locs = find(diff(data(:, 1))==0); % locations of duplicates
        locs = locs + 1; % we will remove latter of duplicates
        data(locs, :) = []; % remove rows with duplicate values

        % Convert raw data to waypoint script
        lla = data(:, [1 4 10 11]); %lat/long/altitude at each time
        % Express time&lat/long data relative to initial values, set to 0:
        dlla = lla;
        dlla(:, 1) = lla(:, 1) - lla(1, 1);
        dlla(:, 3) = lla(:, 3) - lla(1, 3);
        dlla(:, 4) = lla(:, 4) - lla(1, 4);


        while 1 % Add any part of trajectory > desired time
            % End loop if not enough traj remains
            if (dlla(end,1) < minT)
                break;
            end
            % Break up trajectories that are too long
            if (dlla(end,1) > maxT)
                loc = find(dlla(:,1)>=maxT,1);
            else
                loc = length(dlla(:,1));
            end
            dllaT = dlla(1:loc,:); % temporary dlla
            llaT = lla(1:loc,:); % temporary lla

            % Convert lat/long to feet relative to initial position of (0,0)
            pos = dllaT;
            pos(:, 3) = R*dllaT(:, 3)*pi/180; % north/south distance
            pos(:, 4) = R*cosd(.5*(llaT(:,3)+llaT(1,3))).*dllaT(:,4)*pi/180; % east/west distance

            % Scale based on desired speed
            d = diff(pos); % difference at each time step
            r = sqrt(d(:,3).^2+d(:,4).^2); % total distance traveled
            rDes = d(:,1)*speed; % desired total distance traveled
            ratio = rDes ./ r; % ratio
            dadj = [d(:,1:2) d(:,3).*ratio, d(:,4).*ratio]; % adjusted difference
            posAdj = [pos(1,:); pos(1,:)+cumsum(dadj)]; % adjusted position

            % Use data to create waypoints
            waypoints(1, index).initial = posAdj(1, [3 4 2])';
            waypoints(1, index).update = posAdj(2:end, [1 3 4 2])';

            a = posAdj; % makes next few lines shorter
            if numAc == 2 % create dummy straight-line intruder
                waypoints(2, index).initial = [a(end,3) a(end,4) a(1,2)]';
                waypoints(2, index).update = [a(end,1) a(1,3) a(1,4) a(end,2)]';
            end

            % Increment array index if current file was used
            index = index + 1;

            % Remove used portion from trajectory
            if loc == length(lla(:,1)) % last point used was end of traj
                break
            end
            lla = lla(loc:end,:);
            dlla = dlla(loc:end,:);
            dlla(:, 1) = dlla(:, 1) - dlla(1, 1);
            dlla(:, 3) = dlla(:, 3) - dlla(1, 3);
            dlla(:, 4) = dlla(:, 4) - dlla(1, 4);

        end


        %     % Remove trajectories that are too short
        %     if (dlla(end,1) < minT)
        %         fprintf('Skipping encounter %g. Time: %g\n', i,dlla(end,1));
        %         continue; % skip trajectories shorter than desired time
        %     end
        %     % Truncate trajectories that are too long
        %     if (dlla(end,1) > maxT)
        %         loc = find(dlla(:,1)>=maxT,1);
        %         dlla = dlla(1:loc,:);
        %         lla = lla(1:loc,:);
        %     end
        %     % Convert lat/long to feet relative to initial position of (0,0)
        %     pos = dlla;
        %     pos(:, 3) = R*dlla(:, 3)*pi/180; % north/south distance
        %     pos(:, 4) = R*cosd(.5*(lla(:,3)+lla(1,3))).*dlla(:,4)*pi/180; % east/west distance
        %
        %     % Scale based on desired speed
        %     d = diff(pos); % difference at each time step
        %     r = sqrt(d(:,3).^2+d(:,4).^2); % total distance traveled
        %     rDes = d(:,1)*speed; % desired total distance traveled
        %     ratio = rDes ./ r; % ratio
        %     dadj = [d(:,1:2) d(:,3).*ratio, d(:,4).*ratio]; % adjusted difference
        %     posAdj = [pos(1,:); pos(1,:)+cumsum(dadj)]; % adjusted position
        %
        %     % Use data to create waypoints
        %     waypoints(1, index).initial = posAdj(1, [3 4 2])';
        %     waypoints(1, index).update = posAdj(2:end, [1 3 4 2])';
        %
        %     a = posAdj; % makes next few lines much shorter
        %     if numAc == 2 % create dummy straight-line intruder
        %         waypoints(2, index).initial = [a(end,3) a(end,4) a(1,2)]';
        %         waypoints(2, index).update = [a(end,1) a(1,3) a(1,4) a(end,2)]';
        %     end
        %
        %     % Increment array index if current file was used
        %     index = index + 1;

    end

    % Create waypoint and script files
    save_waypoints(foutWP, waypoints);
%     scripts = waypoint2script(waypoints);
%     % ensure every last update has zero accel
%     for i = 1:length(scripts)
%         tf = scripts(1,i).update(1,end);
%         scripts(1,i).update = [scripts(1,i).update [tf+1 0 0 0]'];
%         scripts(1,i).update(3,1) = 0;
%     end
%     save_scripts(foutS, scripts);
end

function pair_uas2_wp(input, fout, numTrials, minSpeedOwn, minSpeedInt, tcpa, speedDistributionFile, runPar)
    % Note that the output file will have both aircraft always starting at the
    % same position. I use additional scripts to apply desired miss distances
    % and perform other postprocessing. 

    % Ths script should support parallel processing if the lines under
    % "Initialize parallelization" are uncommented and the for loops are
    % changed to parfor loops
    % sUAS speed distribution
    data = load(speedDistributionFile); % gives 'cdfVals', 'cdfSpeeds'
    cdfSpeeds = data.cdfSpeeds;
    cdfVals = data.cdfVals;
    
    % sUAS trajectories
    wp = load_waypoints(input);
    wp = wp(1,:);
    numTraj = length(wp);
    
    % Initialize
    ownInitialAll = cell(numTrials,1);
    ownUpdateAll = cell(numTrials,1);
    intInitialAll = cell(numTrials,1);
    intUpdateAll = cell(numTrials,1);
    randInt = randi(numTraj,numTrials,1);
    randOwn = randi(numTraj,numTrials,1);
    
    if runPar
        % Initialize parallelization
        poolobj= gcp('nocreate');
        if(~isempty(poolobj))
            %     delete(gcp('nocreate'));
            %     parpool();
        else
            parpool();
        end
        parfor n = 1:numTrials
            
            % Intruder parameters (randomly selected)
            trajInt = randInt(n);
            intInitial = wp(trajInt).initial;
            intUpdates = wp(trajInt).update;
            
            % Ownship parameters (random from sUAS trajectories)
            trajOwn = randOwn(n);
            ownInitial = wp(trajOwn).initial;
            ownUpdates = wp(trajOwn).update;
            
            % Find sUAS speeds (random)
            while 1
                sInd = find(rand<cdfVals,1); % random index for speed
                if isempty(sInd) % occurs if rand is too close to 1
                    sInd = length(cdfVals);
                end
                % assign speed based on CDF and random index
                ownSpeed = rand*(cdfSpeeds(sInd)-cdfSpeeds(sInd-1)) + cdfSpeeds(sInd-1);
                if ownSpeed > minSpeedOwn
                    break
                end
            end
            ownSpeed = ownSpeed * 1.68781; % kts to ft/s
            % Intruder speed
            while 1
                sInd = find(rand<cdfVals,1); % random index for speed
                if isempty(sInd) % occurs if rand is too close to 1
                    sInd = length(cdfVals);
                end
                % assign speed based on CDF and random index
                intSpeed = rand*(cdfSpeeds(sInd)-cdfSpeeds(sInd-1)) + cdfSpeeds(sInd-1);
                if intSpeed > minSpeedInt
                    break
                end
            end
            intSpeed = intSpeed * 1.68781; % kts to ft/s
            
            % Scale sUAS x-y coords based on desired speed
            pos = [0 ownInitial(1:2)'; ownUpdates([1 2 3],:)'];
            d = diff(pos); % difference at each time step
            r = sqrt(d(:,2).^2+d(:,3).^2); % total distance traveled
            rDes = d(:,1)*ownSpeed; % desired total distance traveled
            ratio = rDes ./ r; % ratio
            dadj = [d(:,1) d(:,2).*ratio, d(:,3).*ratio]; % adjusted difference
            posAdj = [pos(1,:); pos(1,:)+cumsum(dadj)]; % adjusted position
            ownInitial = [posAdj(1,2:3)'; ownInitial(3)];
            ownUpdates = [posAdj(2:end,:)'; ownUpdates(4,:)];
            
            % Repeat for intruder
            pos = [0 intInitial(1:2)'; intUpdates([1 2 3],:)'];
            d = diff(pos); % difference at each time step
            r = sqrt(d(:,2).^2+d(:,3).^2); % total distance traveled
            rDes = d(:,1)*intSpeed; % desired total distance traveled
            ratio = rDes ./ r; % ratio
            dadj = [d(:,1) d(:,2).*ratio, d(:,3).*ratio]; % adjusted difference
            posAdj = [pos(1,:); pos(1,:)+cumsum(dadj)]; % adjusted position
            intInitial = [posAdj(1,2:3)'; intInitial(3)];
            intUpdates = [posAdj(2:end,:)'; intUpdates(4,:)];
            
            % Randomly sample intruder rotation angle
            randAngle = rand*2*pi;
            % Rotate intruder points about origin (where all tracks should start)
            xOld = intUpdates(3,:); % old x (east) points
            yOld = intUpdates(2,:); % old y (north) points
            intUpdates(3,:) = xOld.*cos(randAngle) - yOld.*sin(randAngle);
            intUpdates(2,:) = yOld.*cos(randAngle) + xOld.*sin(randAngle);
            
            % Find intruder position at time of cpa
            index = find(intUpdates(1,:)==tcpa);
            if ~isempty(index) % position update exists at tcpa
                intNcpa = intUpdates(2,index);
                intEcpa = intUpdates(3,index);
                intAcpa = intUpdates(4,index);
            else % interpolate to approximate position
                totalTraj = [[0; intInitial] intUpdates]; % includes t=0 point
                intNcpa = interp1(totalTraj(1,:),totalTraj(2,:),tcpa,'pchip');
                intEcpa = interp1(totalTraj(1,:),totalTraj(3,:),tcpa,'pchip');
                intAcpa = interp1(totalTraj(1,:),totalTraj(4,:),tcpa,'pchip');
            end
            
            % Find o position at time of cpa
            index = find(ownUpdates(1,:)==tcpa);
            if ~isempty(index) % position update exists at tcpa
                ownNcpa = ownUpdates(2,index);
                ownEcpa = ownUpdates(3,index);
                ownAcpa = ownUpdates(4,index);
            else % interpolate
                totalTraj = [[0; ownInitial] ownUpdates]; % includes t=0 point
                ownNcpa = interp1(totalTraj(1,:),totalTraj(2,:),tcpa,'pchip');
                ownEcpa = interp1(totalTraj(1,:),totalTraj(3,:),tcpa,'pchip');
                ownAcpa = interp1(totalTraj(1,:),totalTraj(4,:),tcpa,'pchip');
            end
            
            % Move intruder position and uas alt to line up trajectories at tcpa
            intInitial(1) = intInitial(1) + (ownNcpa-intNcpa);
            intInitial(2) = intInitial(2) + (ownEcpa-intEcpa);
            intUpdates(2,:) = intUpdates(2,:) + (ownNcpa-intNcpa);
            intUpdates(3,:) = intUpdates(3,:) + (ownEcpa-intEcpa);
            ownInitial(3) = ownInitial(3) + (intAcpa-ownAcpa);
            ownUpdates(4,:) = ownUpdates(4,:) + (intAcpa-ownAcpa);
            
            % Output for waypoint variables
            ownInitialAll{n} = ownInitial';
            ownUpdateAll{n} = ownUpdates;
            intInitialAll{n} = intInitial';
            intUpdateAll{n} = intUpdates;
            
        end
    else
        for n = 1:numTrials
            
            % Intruder parameters (randomly selected)
            trajInt = randInt(n);
            intInitial = wp(trajInt).initial;
            intUpdates = wp(trajInt).update;
            
            % Ownship parameters (random from sUAS trajectories)
            trajOwn = randOwn(n);
            ownInitial = wp(trajOwn).initial;
            ownUpdates = wp(trajOwn).update;
            
            % Find sUAS speeds (random)
            while 1
                sInd = find(rand<cdfVals,1); % random index for speed
                if isempty(sInd) % occurs if rand is too close to 1
                    sInd = length(cdfVals);
                end
                % assign speed based on CDF and random index
                ownSpeed = rand*(cdfSpeeds(sInd)-cdfSpeeds(sInd-1)) + cdfSpeeds(sInd-1);
                if ownSpeed > minSpeedOwn
                    break
                end
            end
            ownSpeed = ownSpeed * 1.68781; % kts to ft/s
            % Intruder speed
            while 1
                sInd = find(rand<cdfVals,1); % random index for speed
                if isempty(sInd) % occurs if rand is too close to 1
                    sInd = length(cdfVals);
                end
                % assign speed based on CDF and random index
                intSpeed = rand*(cdfSpeeds(sInd)-cdfSpeeds(sInd-1)) + cdfSpeeds(sInd-1);
                if intSpeed > minSpeedInt
                    break
                end
            end
            intSpeed = intSpeed * 1.68781; % kts to ft/s
            
            % Scale sUAS x-y coords based on desired speed
            pos = [0 ownInitial(1:2)'; ownUpdates([1 2 3],:)'];
            d = diff(pos); % difference at each time step
            r = sqrt(d(:,2).^2+d(:,3).^2); % total distance traveled
            rDes = d(:,1)*ownSpeed; % desired total distance traveled
            ratio = rDes ./ r; % ratio
            dadj = [d(:,1) d(:,2).*ratio, d(:,3).*ratio]; % adjusted difference
            posAdj = [pos(1,:); pos(1,:)+cumsum(dadj)]; % adjusted position
            ownInitial = [posAdj(1,2:3)'; ownInitial(3)];
            ownUpdates = [posAdj(2:end,:)'; ownUpdates(4,:)];
            
            % Repeat for intruder
            pos = [0 intInitial(1:2)'; intUpdates([1 2 3],:)'];
            d = diff(pos); % difference at each time step
            r = sqrt(d(:,2).^2+d(:,3).^2); % total distance traveled
            rDes = d(:,1)*intSpeed; % desired total distance traveled
            ratio = rDes ./ r; % ratio
            dadj = [d(:,1) d(:,2).*ratio, d(:,3).*ratio]; % adjusted difference
            posAdj = [pos(1,:); pos(1,:)+cumsum(dadj)]; % adjusted position
            intInitial = [posAdj(1,2:3)'; intInitial(3)];
            intUpdates = [posAdj(2:end,:)'; intUpdates(4,:)];
            
            % Randomly sample intruder rotation angle
            randAngle = rand*2*pi;
            % Rotate intruder points about origin (where all tracks should start)
            xOld = intUpdates(3,:); % old x (east) points
            yOld = intUpdates(2,:); % old y (north) points
            intUpdates(3,:) = xOld.*cos(randAngle) - yOld.*sin(randAngle);
            intUpdates(2,:) = yOld.*cos(randAngle) + xOld.*sin(randAngle);
            
            % Find intruder position at time of cpa
            index = find(intUpdates(1,:)==tcpa);
            if ~isempty(index) % position update exists at tcpa
                intNcpa = intUpdates(2,index);
                intEcpa = intUpdates(3,index);
                intAcpa = intUpdates(4,index);
            else % interpolate to approximate position
                totalTraj = [[0; intInitial] intUpdates]; % includes t=0 point
                intNcpa = interp1(totalTraj(1,:),totalTraj(2,:),tcpa,'pchip');
                intEcpa = interp1(totalTraj(1,:),totalTraj(3,:),tcpa,'pchip');
                intAcpa = interp1(totalTraj(1,:),totalTraj(4,:),tcpa,'pchip');
            end
            
            % Find o position at time of cpa
            index = find(ownUpdates(1,:)==tcpa);
            if ~isempty(index) % position update exists at tcpa
                ownNcpa = ownUpdates(2,index);
                ownEcpa = ownUpdates(3,index);
                ownAcpa = ownUpdates(4,index);
            else % interpolate
                totalTraj = [[0; ownInitial] ownUpdates]; % includes t=0 point
                ownNcpa = interp1(totalTraj(1,:),totalTraj(2,:),tcpa,'pchip');
                ownEcpa = interp1(totalTraj(1,:),totalTraj(3,:),tcpa,'pchip');
                ownAcpa = interp1(totalTraj(1,:),totalTraj(4,:),tcpa,'pchip');
            end
            
            % Move intruder position and uas alt to line up trajectories at tcpa
            intInitial(1) = intInitial(1) + (ownNcpa-intNcpa);
            intInitial(2) = intInitial(2) + (ownEcpa-intEcpa);
            intUpdates(2,:) = intUpdates(2,:) + (ownNcpa-intNcpa);
            intUpdates(3,:) = intUpdates(3,:) + (ownEcpa-intEcpa);
            ownInitial(3) = ownInitial(3) + (intAcpa-ownAcpa);
            ownUpdates(4,:) = ownUpdates(4,:) + (intAcpa-ownAcpa);
            
            % Output for waypoint variables
            ownInitialAll{n} = ownInitial';
            ownUpdateAll{n} = ownUpdates;
            intInitialAll{n} = intInitial';
            intUpdateAll{n} = intUpdates;
            
        end
    end
    
    waypoints = struct('initial', cell(2, numTrials), 'update', cell(2, numTrials));
    for n = 1:numTrials
        
        waypoints(1,n).initial = ownInitialAll{n};
        waypoints(1,n).update = ownUpdateAll{n};
        waypoints(2,n).initial = intInitialAll{n};
        waypoints(2,n).update = intUpdateAll{n};
        
    end
    
    % Save results
    fprintf('Simulation complete. Beginning save.\n');
    save_waypoints(fout, waypoints);
    fprintf('Save complete\n');
end

function adjWP(input, fout, tcpa, minAlt)
% The reason this is done is because, with scripted encounters, it is hard to know aircraft position in the middle of the encounter without actually simulating
% This script makes the following edits:
% 	Trajectories are moved to have desired miss distances at cpa
%		Vertical distance is sampled from a Gaussian(0,100ft)
%		Horizontal distance is sampled uniformly from +/-5000 ft
%	Trajectories are raised to at least 'minAlt' altitude
%		This ensures TCAS and ACAS X will always issue RAs
%	Trajectories that start extremely close are thrown out
    % Metrics to analyze
    m = cell(1,8);
    names = cell(1,8);
    m{1} = fun.x('x.0',tcpa);
    names{1} = ['x_x_0_' num2str(tcpa)];
    m{2} = fun.x('y.0',tcpa);
    names{2} = ['x_y_0_' num2str(tcpa)];
    m{3} = fun.x('h.0',tcpa);
    names{3} = ['x_h_0_' num2str(tcpa)];
    m{4} = fun.x('x.0',0);
    names{4} = 'x_x_0_0';
    m{5} = fun.x('y.0',0);
    names{5} = 'x_y_0_0';
    m{6} = fun.x('h.0',0);
    names{6} = 'x_h_0_0';
    m{7} = fun.x('x.1',tcpa);
    names{7} = ['x_x_1_' num2str(tcpa)];
    m{8} = fun.x('y.1',tcpa);
    names{8} = ['x_y_1_' num2str(tcpa)];
    m{9} = fun.x('h.1');
    names{9} = 'x_h_1';
    m{10} = fun.x('dx.1',tcpa);
    names{10} = ['x_dx_1_' num2str(tcpa)];
    m{11} = fun.x('dy.1',tcpa);
    names{11} = ['x_dy_1_' num2str(tcpa)];

    c = csim('ss');
    c.encounter = encounter.waypoint_horizontal(input);
    c.encounter.max_steps = 120;

    oldWaypoints = load_waypoints(input);
    fprintf('Loaded waypoints.\n');
    Own = oldWaypoints(1,:);
    Int = oldWaypoints(2,:);
    numEnc = length(Own);
    edits = false(numEnc,1);
    %initX = zeros(numEnc,1);
    %initY = zeros(numEnc,1);
    %initH = zeros(numEnc,1);

    % Offsets
    %horOffset = normrnd(0,0,numEnc,1);
    horOffset = rand(numEnc,1)*5000 - 2500; % unif(-5000, 5000)
    altOffset = normrnd(0,100,numEnc,1);

    % old initial position
    %for i = 1:numEnc
    %    initX(i) = Int(i).initial(3);
    %    initY(i) = Int(i).initial(2);
    %    initH(i) = Int(i).initial(4);
    %end

    % Simulate
    S = c.eval('all', m);

    % Use results to edit encounters
    for n = 1:numEnc
        % Old waypoint values
        ownI = Own(n).initial;
        ownU = Own(n).update;
        intI = Int(n).initial;
        intU = Int(n).update;
        % Find necessary positions
        x0 = S.(names{1})(n);
        y0 = S.(names{2})(n);
        h0 = S.(names{3})(n);
        xi0 = S.(names{4})(n);
        yi0 = S.(names{5})(n);
        hi0 = S.(names{6})(n);
        x1 = S.(names{7})(n);
        y1 = S.(names{8})(n);
        h1 = S.(names{9})(tcpa,n);
        hmin = min(S.(names{9})(:,n));
        dx1 = S.(names{10})(n);
        dy1 = S.(names{11})(n);

        % Adjust (horizontal offset applied tangential to int flight path at cpa)
        ang = atan(dy1/dx1);
        %initY(n) = initY(n) + (y0-y1) + horOffset(n)*cos(ang);
        %initX(n) = initX(n) + (x0-x1) - horOffset(n)*sin(ang);
        %initH(n) = initH(n) + (h0-h1) + altOffset(n);
        intI(1) = intI(1) + (y0-y1) + horOffset(n)*cos(ang);
        intI(2) = intI(2) + (x0-x1) - horOffset(n)*sin(ang);
        intU(2,:) = intU(2,:) + (y0-y1) + horOffset(n)*cos(ang);
        intU(3,:) = intU(3,:) + (x0-x1) - horOffset(n)*sin(ang);
        hmin = hmin + (h0-h1) + altOffset(n);
        if hmin < minAlt % altitude of both aircraft needs to be raised
            %initH(n) = initH(n) + minAlt - hmin;
            %Own(n).initial(4) = Own(n).initial(4) + minAlt - hmin;
            intI(3) = intI(3) + minAlt - hmin;
            intU(4,:) = intU(4,:) + minAlt - hmin;
            ownI(3) = ownI(3) + minAlt - hmin;
            ownU(4,:) = ownU(4,:) + minAlt - hmin;
        end
        % Ensure aircraft do not start too close
        hi1 = intI(3); % new intruder starting height
        yi1 = intI(1); % new intruder starting position
        xi1 = intI(2);
        hi0 = ownI(3); % new ownship starting height
        %minRange = 6070;
        %minVert = 250;
        minSlant = 5000;
        %initRange = sqrt((xi0-xi1)^2+(yi0-yi1)^2);
        initSlant = sqrt((xi0-xi1)^2+(yi0-yi1)^2+(hi0-hi1)^2);
        %if initRange < minRange && abs(hi0-hi1) < minVert
        %    edits(n) = true; % this will cause encounter to be removed
        %end
        if initSlant < minSlant
            edits(n) = true;
        end

        % Recreate waypoints
        Own(n).initial = ownI;
        Own(n).update = ownU;
        Int(n).initial = intI;
        Int(n).update = intU;

    end


    % Recreate intruder script
    %for i = 1:numEnc
    %    Int(i).initial(3) = initX(i);
    %    Int(i).initial(2) = initY(i);
    %    Int(i).initial(4) = initH(i);
    %end

    % Save Results
    waypoints = [Own; Int];
    waypoints = waypoints(:, ~edits);
    save_waypoints(fout, waypoints);
    fprintf('Save complete.\n');
end

function finWP(input, fout, numTrials, tcpa, tol)
    steps = 120; % time steps

    c = csim('ss');
    c.encounter = encounter.waypoint_horizontal(input);
    c.encounter.max_steps = steps;

    t = c.eval('all', fun.tcaMS());
    t = t.tcaMS;

    index = (t<tcpa-tol) + (t>tcpa+tol);

    waypoints = load_waypoints(input);
    waypoints = waypoints(:, ~index);
    numEnc = length(waypoints);
    if numTrials < numEnc % we can't save more encounters than we have
        waypoints = waypoints(:, 1:numTrials); 
    end
    save_waypoints(fout, waypoints);
end
