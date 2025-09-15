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

INTERVAL = 1; % unit: degrees
lng = (40:INTERVAL:64)';
lat = (24:INTERVAL:43)';
[lat,lng] = meshgrid(lat,lng);

wgs84 = wgs84Ellipsoid("meter");

eps_xx = zeros(size(lat)); eps_yy = eps_xx; eps_xy = eps_xx; theta = eps_xx;
omega = eps_xx; dx = eps_xx; dy = eps_xx; eps1 = eps_xx; eps2 = eps_xx; D = eps_xx;
for i = 1:size(lat,2)
    for j = 1:size(lng,1)
        i,j

    P = [lat(1,i) lng(j,1)];
    
    [Dist,Az] = distance(P,[data.("Lat(°N)") data.("Long(°E)")],wgs84);

    % Maximum distance for GPS stations to be used in LS
    THRESHOLD = mean(Dist) + 2.5*std(Dist); % median(Dist)
    index = Dist < THRESHOLD;
    data_used = data(index,:);  Dist = Dist(index);  Az = Az(index);

    Dx = Dist.*sind(Az);
    Dy = Dist.*cosd(Az);
    m = size(Dx,1)*2; n = 6;

    % Design Matrix
    A = zeros(m,n); 
    A(1:2:end,1) = Dx;
    A(2:2:end,2) = Dy;
    A(1:2:end,3) = Dy;
    A(2:2:end,3) = Dx;
    A(1:2:end,4) = Dy;
    A(2:2:end,4) = -Dx;
    A(1:2:end,5) = 1;
    A(2:2:end,6) = 1;

    % Observation Vector
    Y = [data_used.Evel'/1000;data_used.Nvel'/1000];
    Y = reshape(Y,[],1);
    
    % Clustering Weight
    % vlat = [P(1);data_used.("Lat(°N)")]; 
    % vlng = [P(2);data_used.("Long(°E)")];
    % [V,C,ll] = VoronoiLimit(vlng,vlat,"figure","off");
    % vor_area = zeros(m/2+1,1);
    % for k = 1:length(C)
    %     v1 = V(C{k},1) ;
    %     v2 = V(C{k},2) ;
    %     vor_area(k,1) = polyarea(v1,v2);
    % end
    % total_area = sum(vor_area);
    % vor_area(1) = [];
    % Z = size(vor_area,1).*vor_area/total_area;
    % Z = [Z';Z'];
    % Z = reshape(Z,[],1);

    % Distance Weight
    R = Dist/1000; % unit: km
    % for DD = 0:50:5000 % unit: km
    %     L = exp(-R.^2/DD^2);
    %     if sum(L > 0.5) > size(R,1)/2
    %         D(j,i) = DD;
    %         break
    %     end
    % end
    L = exp(-R.^2/200^2);
    L = [L';L'];
    L = reshape(L,[],1);

    % Variance & Covariance Weight
    C = zeros(m);
    C(1:2:end,1:2:end) = diag((data_used.SigVe/1000).^2);
    C(2:2:end,2:2:end) = diag((data_used.SigVn/1000).^2);
    C(2:2:end,1:2:end) = diag((data_used.SigVe/1000).*(data_used.SigVn/1000).*data_used.Cor); 
    C(1:2:end,2:2:end) = diag((data_used.SigVe/1000).*(data_used.SigVn/1000).*data_used.Cor);

    % Total Weight 
    Z = ones(m,1);
    G = L.*Z;
    W = inv(C).*G;

    % Least Squares 
    x_hat = lscov(A,Y,W,"chol");
    eps_xx(j,i) = x_hat(1);
    eps_yy(j,i) = x_hat(2);
    eps_xy(j,i) = x_hat(3);
    omega(j,i) = x_hat(4);
    dx(j,i) = x_hat(5);
    dy(j,i) = x_hat(6);
    
    EPS = [eps_xx(j,i) eps_xy(j,i); eps_xy(j,i) eps_yy(j,i)];
    pEPS = eig(EPS);
    eps1(j,i) = pEPS(2);
    eps2(j,i) = pEPS(1);
    theta(j,i) = 0.5*atan2(2*eps_xy(j,i),(eps_xx(j,i)-eps_yy(j,i)));

    end
end

I1 = eps_xx + eps_yy;
I2 = eps_xx.*eps_yy - eps_xy.^2;
result.Latitude = lat;
result.Longitude = lng;
result.eps_xx = eps_xx;
result.eps_yy = eps_yy;
result.eps_xy = eps_xy;
result.eps1 = eps1;
result.eps2 = eps2;
result.theta = theta;
result.I1 = I1;
result.I2 = I2;
result.omega = omega;
result.dx = dx;
result.dy = dy;
result.D = D;

clearvars -except result data grid

