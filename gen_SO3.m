function [R, alpha, beta, gamma] = gen_SO3(alpha, beta, gamma)
%GEN_SO3    randomly generate a matrix R in SO3
%   g = GEN_SO3
%   g = GEN_SO3(alpha, beta, gamma)

%-- Auther: hshi17 11/17/18 --%

    if nargin == 0
        alpha = 2*pi*rand;
        beta = pi*rand;
        gamma = 2*pi*rand;
    end  
    
    R = eul2rotm([alpha, beta, gamma], 'ZYZ');
    
end