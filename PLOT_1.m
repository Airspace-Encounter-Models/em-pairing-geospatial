%% Inputs
% Parameters from RUN_1.m
iso_3166_2 = {'US-CA','US-KS','US-MA','US-MS','US-NC','US-NV','US-NY','US-TX','US-VA'};
usecase = 'all-noshield';
anchorRange_nm = 1.25; % Distance pairs need to be away from anchor point 30*(150 / 3600) = 1.25...closing speed of 150 knots for 30 seconds

% Edges to discretize anchors.num_geospatial
edges = [0,1,2,4,6,inf];

% Color order
colorOrder = [204 121 167;... % Reddish purple
    230 159 0;... % Orange
    0 114 178;... % Blue
    86 180 233;... % Sky Blue
    0 158 155;... % Bluish green
    213 94 0;.... % Vermillion
    240 228 66;... % Yellow
    ] / 255;

% Map bounds buffer
buff_deg = 0.05;

%% Iterate over locations
for j=1:1:numel(iso_3166_2)
    % Output from RUN_1.m
    jhash = [usecase '_' iso_3166_2{j} '-anchorRange' num2str(anchorRange_nm)];
    inFile = ['output' filesep 'pairs-' jhash '.mat'];
    
    % Load
    load(inFile,'anchors');
    
    % Discretize
    Y = discretize(anchors.num_geospatial,edges);
    
    % Create GeographicAxes and Base Map
    figure(j); set(gcf,'name',inFile);
    gx = geoaxes('Basemap','streets-dark',...
        'FontSize',16,...
        'FontWeight','bold',...
        'Position',[0 0 1 1]);
    bmap = geobasemap(gx);
    
    % Calculate bounds
    %minLat_deg = min(anchors.lat_deg) - buff_deg;
    %maxLat_deg = max(anchors.lat_deg) + buff_deg;
    %minLon_deg = min(anchors.lon_deg) - buff_deg;
    %maxLon_deg = max(anchors.lon_deg) + buff_deg;
    %geolimits(gx,[minLat_deg maxLat_deg],[minLon_deg maxLon_deg]);
    
    % Iterate through discretized anchors.num_geospatial
    hold on
    for i=1:1:max(Y)
        l = Y == i;
        geoscatter(anchors.lat_deg(l), anchors.lon_deg(l),25,colorOrder(i,:),'o','filled');
        
        %geoshow(anchors.lat_deg(l), anchors.lon_deg(l));
        labels{i} = sprintf('[%i, %i)',edges(i),edges(i+1));
    end
    hold off
    legend(labels);
    
    % Adjust map size
    set(gcf,'Units','inches','Position',[1 1 11 8])
    
    % Save figure
    %saveas(gcf,['output' filesep 'plot_1_' jhash '.png'],'png');
    print(['output' filesep 'plot_1_' jhash '.png'],'-dpng','-r300');
end