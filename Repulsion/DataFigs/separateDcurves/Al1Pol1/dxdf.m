function [dxdf, xx] = dxdf(x, y, h)
    % ������ �� x 
    yy = y(1):h:y(length(y));
    % ���������� ����� y
    xx = spline(y,x,yy);
    % ����� ��� ��� � ������ �����
    dxdf = diff(xx)/h;
    n = length(dxdf);
    xx = xx(1:n);

end