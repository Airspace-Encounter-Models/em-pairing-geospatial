function fig = plotEncounter(conflicts,mdlType1,mdlType2,track1,track2,sidx1,eidx1,sidx2,eidx2,lat0_deg,lon0_deg,randSeed)
% Copyright 2019 - 2021, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
% SEE ALSO createEncounters_2

%% Filter based on time
% Filter track timetables
if sidx1 <= eidx1
    track1 = track1(sidx1:1:eidx1,:);
else
    track1 = track1(sidx1:-1:eidx1,:);
end
track1.Properties.StartTime = seconds(0);
track1.Properties.TimeStep = seconds(1);

% Filter track timetables
if sidx2 <= eidx2
    track2 = track2(sidx2:1:eidx2,:);
else
    track2 = track2(sidx2:-1:eidx2,:);
end
track2.Properties.StartTime = seconds(0);
track2.Properties.TimeStep = seconds(1);

%% Parse and calculate
hmd_ft = conflicts(:,9);
vmd_ft = conflicts(:,10);
r_ft = hypot(track1.east_ft-track2.east_ft,track1.north_ft-track2.north_ft);
z_ft = abs(track1.alt_ft_msl - track2.alt_ft_msl);
d_ft = sqrt((track1.north_ft - track2.north_ft).^2 + (track1.east_ft - track2.east_ft).^2 + (track1.alt_ft_msl - track2.alt_ft_msl).^2);

[~,idxCPA] = min(d_ft);

%% Plot
% Create figure
fig = figure('name',sprintf('Encounter seed = %i',randSeed),'Units','inches','Position',[0.5 0.5 16 8]);
tl = tiledlayout(2,5,'Padding','compact');

% Altitude and speed
nexttile(1,[2 1]);
TTdyn = timetable(track1.alt_ft_agl,track1.alt_ft_msl,track2.alt_ft_agl,track2.alt_ft_msl,track1.speed_ft_s,track2.speed_ft_s,track1.speed_ft_s*0.592484,track2.speed_ft_s*0.592484,'RowTimes',track1.Time,...
    'VariableNames',{'AC 1 (ft AGL)','AC 1 (ft MSL)','AC 2 (ft AGL)','AC 2 (ft MSL)','AC 1 (fps)','AC 2 (fps)','AC 1 (kts)','AC 2 (kts)'});
newYlabels = {'Altitude (ft AGL)','Altitude (ft MSL)','Speed (fps)','Speed (knots)'};
stackedplot(TTdyn, {{'AC 1 (ft AGL)','AC 2 (ft AGL)'},{'AC 1 (ft MSL)','AC 2 (ft MSL)'},{'AC 1 (fps)','AC 2 (fps)'},{'AC 1 (kts)','AC 2 (kts)'}},'DisplayLabels',newYlabels)
grid on;

% Distance
nexttile(2,[2 1]);
TTd = timetable(r_ft,z_ft,d_ft,'RowTimes',track1.Time);
stackedplot(TTd,{{'d_ft'},{'r_ft'},{'z_ft'}},'DisplayLabels',{'Distance (ft)','HMD (ft)','VMD (ft)'});
grid on;

% Big plot
nexttile(3,[2 2]);
geoplot(track1.lat_deg,track1.lon_deg,'b-',track2.lat_deg,track2.lon_deg,'r--',...
    track1.lat_deg(1),track1.lon_deg(1),'bs',track2.lat_deg(1),track2.lon_deg(1),'rs');
set(gca,'Basemap','topographic','TickLabelFormat','dd');
hold on;

%Threshold Points
for w=1:1:size(conflicts,1)
    if w==1; hv = 'on'; else hv = 'off'; end; % Only want to plot conflict points once on legend
    geoscatter(conflicts(w,1),conflicts(w,2),60,'b','h','filled','HandleVisibility',hv)
    geoscatter(conflicts(w,3),conflicts(w,4),60,'r','p','filled','HandleVisibility',hv)
    
    % Alternate horizontal alignment for odd / even values of w
    % o1, o2 = offsets for text
    if mod(w,2) == 0
        ha1 = 'left';
        ha2 = 'right';
        o1 = nm2deg(0.04);
        o2 = -nm2deg(0.04);
    else
        ha1 = 'right';
        ha2 = 'left';
        o1 = -nm2deg(0.04);
        o2 = nm2deg(0.04);
    end
    
    text(conflicts(w,1),conflicts(w,2)+o1,num2str(w),'Color','b','FontSize',8,'FontWeight','bold','HorizontalAlignment',ha1);
    text(conflicts(w,3),conflicts(w,4)+o2,num2str(w),'Color','r','FontSize',8,'FontWeight','bold','HorizontalAlignment',ha2);
end
geoplot(lat0_deg,lon0_deg,'k*');
hold off;
lg = legend(sprintf('Aircraft 1 (%s)',mdlType1), sprintf('Aircraft 2 (%s)',mdlType2),'AC1 @ initial','AC2 @ initial','AC1 @ conflict','AC2 @ conflict',sprintf('(%0.3f,%0.3f)',lat0_deg,lon0_deg),...
    'NumColumns',2,'Location','best');

% Local Cartessian (north / east) plot
nexttile(5,[1 1]);
plot(track1.east_ft,track1.north_ft,'b-',track2.east_ft,track2.north_ft,'r--',0,0,'k*');
grid on; axis square; axis equal; xlabel('East (ft)'); ylabel('North (ft)');
%legend(sprintf('Aircraft 1 (%s)',mdlType1), sprintf('Aircraft 2 (%s)',mdlType2),'(0,0)');

% nexttile(8,[1 1]);
% geoplot(track1.lat_deg,track1.lon_deg,'b-',track2.lat_deg,track2.lon_deg,'r--',lat0_deg,lon0_deg,'k*');
% set(gca,'TickLabelFormat','dd');
%legend(sprintf('Aircraft 1 (%s)',mdlType1), sprintf('Aircraft 2 (%s)',mdlType2),sprintf('(%0.3f,%0.3f)',lat0_deg,lon0_deg));

% HMD / VMD
nexttile(10,[1 1]);
scatter(r_ft(1),z_ft(1),36*2, 'g','s','filled','DisplayName','Initial'); hold on;
scatter(r_ft(idxCPA),z_ft(idxCPA),36*2, 'r','o','filled','DisplayName','CPA'); 
scatter(r_ft(end),z_ft(end),36*2, 'b','d','filled','DisplayName','Final'); 
scatter(hmd_ft,vmd_ft,36*2, 'k','p','filled','DisplayName','Conflict');hold off;

legend('Location','best');
axis([0 round(max(r_ft),-2)+100 0 round(max(z_ft),-2)+100]);
grid on; xlabel('HMD (ft)'); ylabel('VMD (ft)');

% Legend of big plot
%lg.Layout.Tile = 'south'; % <-- place legend south of tiles
%lg.NumColumns = 4;

% Adjust map size, legend, hold off
% print(['output' filesep 'plot_encounter_examples' filesep 'example_2' '_lk' num2str(lk(i)) '_k' num2str(k) '_i' num2str(i) '_lat' num2str(lat0(lk(i))) '_lon' num2str(lon0(lk(i))) '_thresHorz' num2str(thresHorz_ft) '_thresVert' num2str(thresVert_ft) '_ac1' mdlType1 '_ac2' mdlType2 '.png'],'-dpng','-r300');
