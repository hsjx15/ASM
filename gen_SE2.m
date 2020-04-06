function [g, r, phi, theta] = gen_SE2(r, phi, theta)
%GEN_SE2    randomly generate a matrix g in SE2
%   g = GEN_SE2
%   g = GEN_SE2(r)
%   g = GEN_SE2(r, phi, theta)

%-- Auther: hshi17 11/17/18 --%

    if nargin <= 1
        p = 0;
        if nargin == 1
            p = r;
        end
        r = randi(10^p) + rand - 1;  % between 0 to 10^p
        phi = 2*pi*rand;
        theta = 2*pi*rand;
    end  
    
    g = [cos(theta), -sin(theta), r*cos(phi);
         sin(theta), cos(theta), r*sin(phi);
         0, 0, 1];

end