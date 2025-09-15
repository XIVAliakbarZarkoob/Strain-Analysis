%% Aliakbar Zarkoob; AKA "XIV" 810302065

clc, clear, close all, format long g, beep off

%% Load data & Initialization

data = readtable("GPS_data_GeoDynamics.xlsx","NumHeaderLines",2,"FileType","spreadsheet", ...
    "VariableNamingRule","preserve","TextType","string");

data.("Long(°E)") = data.("Long(°E)")/10000;
data.("Lat(°N)") = data.("Lat(°N)")/10000;
% Velocity Units: mm/year
data.("Evel") = data.Evel.double;
data.("Nvel") = data.Nvel.double;
data.("SigVe") = data.SigVe.double;
data.("SigVn") = data.SigVn.double;
data.("Cor") = data.Cor.double;

%% Main

INTERVAL = 0.5; % unit: degrees
LAT = (24:INTERVAL:43)';
LNG = (40:INTERVAL:64)';
[LAT,LNG] = meshgrid(LAT,LNG);

lat = LAT(:); lng = LNG(:);
gps_num = height(data);
m = gps_num*2; % total observation number
n = 6;

% GPS stations Latitude and Longitude 
% (used for distance and azimuth calculations)
dlat = repmat(data.("Lat(°N)"),size(lat,1),1);
dlng = repmat(data.("Long(°E)"),size(lat,1),1);

% Grid Latitude and Longitude 
% (used for distance and azimuth calculations)
lat = repelem(lat,gps_num);
lng = repelem(lng,gps_num);

% Calculating Distance and Azimuth from each point in grid to each GPS
% station (with wgs84 as the reference ellipsoid)
wgs84 = wgs84Ellipsoid("meter");
[Dist,Az] = distance([lat lng],[dlat dlng],wgs84);
Dx = Dist.*sind(Az);
Dy = Dist.*cosd(Az);

% Variance-Covariance matrix of observations
C = zeros(m);
C(1:2:end,1:2:end) = diag((data.SigVe/1000).^2);
C(2:2:end,2:2:end) = diag((data.SigVn/1000).^2);
C(2:2:end,1:2:end) = diag((data.SigVe/1000).*(data.SigVn/1000).*data.Cor);
C(1:2:end,2:2:end) = diag((data.SigVe/1000).*(data.SigVn/1000).*data.Cor);

% Design matrix
A = zeros(size(Dist,1)*2,n);
A(1:2:end,1) = Dx;
A(2:2:end,2) = Dy;
A(1:2:end,3) = Dy;
A(2:2:end,3) = Dx;
A(1:2:end,4) = Dy;
A(2:2:end,4) = -Dx;
A(1:2:end,5) = 1;
A(2:2:end,6) = 1;

% Observation vector
Y = [data.Evel';data.Nvel']/1000;
Y = reshape(Y,[],1);

R = repelem(Dist/1000,2);


vlat = data.("Lat(°N)");

vlng = data.("Long(°E)");

[V,P,ll] = VoronoiLimit(vlng,vlat,"figure","off");
vor_area = zeros(size(P,1),1);
for k = 1:length(P)
    v1 = V(P{k},1) ;
    v2 = V(P{k},2) ;
    vor_area(k,1) = polyarea(v1,v2);
end
Z = size(vor_area,1).*vor_area/sum(vor_area);
Z = repelem(Z,2);
x_hat = zeros(6,size(LAT(:),1));

for i = 1:size(LAT(:),1)


    DD = repelem(Dist(i*gps_num-(gps_num-1):i*gps_num),2);
    THRESHOLD = mean(DD) + 2.5*std(DD);
    index = DD > THRESHOLD;

    ZZ = Z;
    ZZ(index) = [];

    D =  (mean(DD) + 1*std(DD))/10000;
    L = exp(-R.^2/D^2);

    AA = A(i*m-(m-1):i*m,:);
    AA(index,:) = [];
    LL = L(i*m-(m-1):i*m);
    LL(index) = [];

    CC = C;
    CC(index,:) = [];
    CC(:,index) = [];
    W = inv(CC).*LL.*ZZ;
    YY = Y;
    YY(index) = [];

    x_hat(:,i) = inv(AA'*W*AA)*AA'*W*YY;

end

eps_xx = x_hat(1,:);
eps_yy = x_hat(2,:);
eps_xy = x_hat(3,:);

figure
hold on
worldmap([24,43],[40,64])
coast = load('coastlines.mat');
contourfm(LAT,LNG,reshape(eps_xx,size(LAT))*10^9)
plotm(coast.coastlat,coast.coastlon,"k",LineWidth=2)

% geoshow('landareas.shp', 'FaceColor', [1 1 1])







