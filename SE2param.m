function [r, phi, theta] = SE2param(g)
%SE22PARAM
%
%   Input:
%    - g: 3*3 matrix in SE(2)
%   Output:
%    - r: a real number
%    - phi: between 0 and 2*pi
%    - theta: between 0 and 2*pi

%-- Auther: hshi17 11/17/18 --%

    r = sqrt(g(1,3)^2 + g(2,3)^2);
    
    phi = atan2(g(2,3), g(1,3));
    if phi < 0
        phi = phi + 2*pi;   % between 0 and 2*pi
    end
    
    theta = atan2(g(2,1), g(1,1));
    if theta < 0
        theta = theta + 2*pi;   % between 0 and 2*pi
    end
    
end