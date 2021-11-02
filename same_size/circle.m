function h = circle(x,y,r)
% draws a circle with radius r around x and y
th = 0:pi/50:2*pi;
x_unit = r * cos(th) + x;
y_unit = r * sin(th) + y;
%h = plot(xunit, yunit,'k');
h=fill(x_unit, y_unit, 'k');
