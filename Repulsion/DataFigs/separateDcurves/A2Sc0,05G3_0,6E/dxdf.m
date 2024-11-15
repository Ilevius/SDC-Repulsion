function [dxdf, xx] = dxdf(x, y, h)
    % шагаем по x 
    yy = y(1):h:y(length(y));
    % интерполир соотв y
    xx = spline(y,x,yy);
    % произ обр фун в каждой точке
    dxdf = diff(xx)/h;
    n = length(dxdf);
    xx = xx(1:n);

end