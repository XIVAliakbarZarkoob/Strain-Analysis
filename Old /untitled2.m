clc
clear
close all

rng default;
x = rand([1 10]);
y = rand([1 10]);

voronoi(x,y)

[V,C,XY] = VoronoiLimit(x',y');
% hold on
% plot(XY,'*')

for i = 1:length(C)
    v1 = V(C{i},1) ;
    v2 = V(C{i},2) ;
    idx = find(inpolygon(XY(:,1),XY(:,2),v1,v2));
    % hold on
    pgon = polyshape(v1,v2);
    % plot(pgon)
    vorarea(i,1:2) = [polyarea(v1,v2) ] ;
end