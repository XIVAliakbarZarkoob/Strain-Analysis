% XIV & XERO
function ellipse(xc,yc,a,b,theta,color,linewidth)

t=-10:0.01:10;
X=a*cos(t); X=X';
Y=b*sin(t); Y=Y';
R=[cos(theta) -sin(theta);sin(theta) cos(theta)];

for i=1:length(t)

    XY=R*[X(i);Y(i)];
    X(i)=XY(1); Y(i)=XY(2);

end
X=X+xc;
Y=Y+yc;
plotm(X,Y,"Color",color,"LineWidth",linewidth)


end