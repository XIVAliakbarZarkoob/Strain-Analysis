%% Plot & Maps

geobasemap darkwater
geoplot(data.("Lat(°N)"),data.("Long(°E)"),"k^","MarkerFaceColor","k","MarkerSize",5)

%% GPS velocity field
figure
hold on
worldmap([22,45],[38,66])
geoshow('landareas.shp', 'FaceColor', [1 1 1])
irboarder = readtable('ir.csv'); 
irboarder.Properties.VariableNames =["Latitude" "Longitude"];
plotm(irboarder.Latitude, irboarder.Longitude, 'k.',"MarkerSize",4)
quiverm(data.("Lat(°N)"), data.("Long(°E)"), data.Nvel, data.Evel,'b',"off")
title('GPS velocity field relative to the Eurasia fixed frame')

%% 
figure
hold on
worldmap([22,45],[38,66])
geoshow('landareas.shp', 'FaceColor', [1 1 1])
contourfm(result.Latitude,result.Longitude,result.eps_xx)
plotm(irboarder.Latitude, irboarder.Longitude, 'k.',"MarkerSize",4)
hcb = colorbar;
set(get(hcb,'Ylabel'),'String','Exx strain in nanostrain / year')
colormap jet

%%
worldmap([22,45],[38,66])
geoshow('landareas.shp', 'FaceColor', [1 1 1])
surfm(grid(:,1),grid(:,2),result.eps_xx)
colormap("jet")


%% Strain Ellipse

worldmap([22,45],[38,66])
geoshow('landareas.shp', 'FaceColor', [1 1 1])
hold on
for i = 1:size(grid,1)

    SCALE = 5*10^7;
    c = grid(i,:)';
    eps1 = result.eps1(i)*SCALE;
    eps2 = result.eps2(i)*SCALE;
    theta = result.theta(i);

    ax1s = [-abs(eps1);0]; ax1e = [abs(eps1);0];
    ax2s = [0;-abs(eps2)]; ax2e = [0;abs(eps2)];
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];

    ax1s = R*ax1s; ax1e = R*ax1e;
    ax2s = R*ax2s; ax2e = R*ax2e;
    ax1s = c+ax1s; ax1e = c+ax1e;
    ax2s = c+ax2s; ax2e = c+ax2e;

    if eps1 > 0 
        color1 = "r";
    else
        color1 = "b";
    end
    if eps2 > 0
        color2 = "r";
    else
        color2 = "b";
    end
    % Red: Tension  &  Blue: Compression
    plotm([ax1s(1) ax1e(1)],[ax1s(2) ax1e(2)],"Color",color1,"LineWidth",2)
    plotm([ax2s(1) ax2e(1)],[ax2s(2) ax2e(2)],"Color",color2,"LineWidth",2)



end


