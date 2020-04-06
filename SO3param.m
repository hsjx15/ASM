function [alpha, beta, gamma] = SO3param(R)
%SO32PARAM
%
%   Input:
%    - g: 3*3 matrix in SO(3)
%   Output:
%    - alpha: between 0 and 2*pi
%    - beta: between 0 and pi
%    - gamma: between 0 and 2*pi

%-- Auther: hshi17 11/17/18 --%

    beta = acos(R(3,3));    % between 0 and pi

    if R(3,3) ~= 1
        alpha = atan2(R(2,3), R(1,3));
        if alpha < 0
            alpha = alpha + 2*pi;   % between 0 and 2*pi
        end

        gamma = atan2(R(3,2), -R(3,1));
        if gamma < 0
            gamma = gamma + 2*pi;   % between 0 and 2*pi
        end
    else
        alpha = atan2(R(2,1), R(1,1))/2;
        gamma = alpha;
    end

end