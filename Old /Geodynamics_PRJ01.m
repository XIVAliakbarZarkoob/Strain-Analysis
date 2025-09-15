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
lng = (40:INTERVAL:64)'; lngn = size(lng,1);
lat = (24:INTERVAL:43)'; latn = size(lat,1);
grid = zeros(latn,2*lngn);
grid(:,1:lngn) = lat.*ones(latn,lngn);
grid(:,lngn+1:end) = (lng.*ones(lngn,latn))';
grid = reshape(grid,[],2); % Latitude Longitude

wgs84 = wgs84Ellipsoid("meter");

eps_xx = zeros(size(grid,1),1); eps_yy = eps_xx; eps_xy = eps_xx; theta = eps_xx;
omega = eps_xx; dx = eps_xx; dy = eps_xx; eps1 = eps_xx; eps2 = eps_xx; D = eps_xx;
for i = 1:size(grid,1)

    P = grid(i,:);
    
    [Dist,Az] = distance(P,[data.("Lat(°N)") data.("Long(°E)")],wgs84);

    % Maximum distance for GPS stations to be used in LS
    THRESHOLD = median(Dist); 
    index = Dist < THRESHOLD;
    data_used = data(index,:);  Dist = Dist(index);  Az = Az(index);

    Dx = Dist.*sind(Az);
    Dy = Dist.*cosd(Az);
    m = size(Dx,1)*2; n = 6;

    % Design Matrix
    A = zeros(m,n); AA = A;
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
    vlat = [P(1);data_used.("Lat(°N)")]; 
    vlng = [P(2);data_used.("Long(°E)")];
    [V,C,ll] = VoronoiLimit(vlng,vlat,"figure","off");
    vor_area = zeros(m/2+1,1);
    for j = 1:length(C)
        v1 = V(C{j},1) ;
        v2 = V(C{j},2) ;
        vor_area(j,1) = polyarea(v1,v2);
    end
    total_area = sum(vor_area);
    vor_area(1) = [];
    Z = size(vor_area,1).*vor_area/total_area;
    Z = [Z';Z'];
    Z = reshape(Z,[],1);

    % Distance Weight
    R = sqrt(Dx.^2 + Dy.^2)/1000; % unit: km
    for DD = 0:50:5000 % unit: km
        L = exp(-R.^2/DD^2);
        if sum(L > 0.5) > size(R,1)/2
            D(i) = DD;
            break
        end
    end
    L = [L';L'];
    L = reshape(L,[],1);

    % Variance & Covariance Weight
    C = zeros(m);
    C(1:2:end,1:2:end) = diag((data_used.SigVe/1000).^2);
    C(2:2:end,2:2:end) = diag((data_used.SigVn/1000).^2);
    C(2:2:end,1:2:end) = diag((data_used.SigVe/1000).*(data_used.SigVn/1000).*data_used.Cor); 
    C(1:2:end,2:2:end) = diag((data_used.SigVe/1000).*(data_used.SigVn/1000).*data_used.Cor);

    % Total Weight 
    G = L.*Z;
    W = inv(C).*G;

    % Least Squares 
    x_hat = lscov(A,Y,W,"chol");
    eps_xx(i) = x_hat(1);
    eps_yy(i) = x_hat(2);
    eps_xy(i) = x_hat(3);
    omega(i) = x_hat(4);
    dx(i) = x_hat(5);
    dy(i) = x_hat(6);
    
    EPS = [eps_xx(i) eps_xy(i); eps_xy(i) eps_yy(i)];
    pEPS = eig(EPS);
    eps1(i) = pEPS(2);
    eps2(i) = pEPS(1);
    theta(i) = 0.5*atan(2*eps_xy(i)/(eps_xx(i)-eps_yy(i)));

end

I1 = eps_xx + eps_yy;
I2 = eps_xx.*eps_yy - eps_xy.^2;
result = table(grid(:,1),grid(:,2),eps_xx,eps_yy,eps_xy,eps1,eps2,theta,I1,I2,omega,dx,dy,D);
result.Properties.VariableNames = ["Latitude","Longitude","eps_xx","eps_yy","eps_xy","eps1",...
    "eps2","theta","I1","I2","omega","dx","dy","D"];

clearvars -except result data grid

