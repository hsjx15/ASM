function u = diff_iur_se3(X, p, s, l1, l, m1, m)
%DIFF_IUR_SE3   Differentiation of IUR for SE(3)
%   u = diff_iur_se3(X, p, s, l1, l);
%   u = diff_iur_se3(X, p, s, l1, l, m1, m);
%
%   Input:
%   - X: 4x4 matrix in se(3)
%   - p: positive real number, first dual index
%   - s: integre, second dual index
%    - l1, m1, l, m: integers, indices of U, where
%           l1, l \ge abs(s), 
%           m1 \in {-l1, l1}, 
%           m \in {-l, l}
%   Output:
%    - u: differentiation of IUR for SE(3), 11*m1*l2*m2 matrix

%-- Auther: hshi17 12/20/18 --%

    if l1(1) < abs(s)
        disp('value of l1 is incorrect');
        u = NaN; return;
    end
    if l(1) < abs(s)
        disp('value of l is incorrect');
        u = NaN; return;
    end
    
    if nargin ~= 5 && nargin ~= 7
        disp('input error');
        u = NaN; return;
    end

    x = se32vec(X);
    
    L = max(l1(end), l(end));
    u_size = (L+1)^2 - s^2;
    u = zeros(u_size, u_size);
    
    if nargin == 7  % l1, l, m1, m are index numbers
        if m1 > l1
            disp('value of m1 is incorrect');
            u = 0; return;
        end
        if m > l
            disp('value of m is incorrect');
            u = 0; return;
        end
        u = diff_iur_se3_one(x, p, s, l1, l, m1, m);
    end
    
    for i = 1:u_size
        for j = 1:u_size
            lc = ceil(sqrt([i j]+s^2))-1;   % current l1 and l
            
            if (lc(1) <= l1(end) && lc(2) <= l(end))    % in boundary
                mc = [i j] - lc .* (lc+1) -1 + s^2;     % current m1 and m
                
                u(i,j) = ...
           diff_iur_se3_one(x, p, s, lc(1), lc(2), mc(1), mc(2));
            end
        end
    end
    
end

function u = diff_iur_se3_one(x, p, s, l1, l, m1, m)

    if l == 0
        u = 0;
    else

        c_mN = c_func(l, -m);
        c_mP = c_func(l, m);

        sigma_l     =   (l1==l);
        sigma_l1N   =   ((l1-1)==l);
        sigma_l1P   =   ((l1+1)==l);
    %     sigma_lN    =   (l1==(l-1));
    %     sigma_lP    =   (l1==(l+1));
        sigma_m     =   (m1==m);
        sigma_m1N   =   ((m1-1)==m);
        sigma_m1P   =   ((m1+1)==m);
        sigma_mN    =   (m1==(m-1));
        sigma_mP    =   (m1==(m+1));

        u1 = -1i/2 * c_mN .* sigma_l .* sigma_m1P ...
            - 1i/2 * c_mP .* sigma_l .* sigma_m1N;
        u2 = +1/2 * c_mN .* sigma_l .* sigma_m1P ...
            - 1/2 * c_mP .* sigma_l .* sigma_m1N;
        u3 = -1i * m .* sigma_l .* sigma_m;

        gamma_m1N = gamma_func(s, l1, -m1);
        gamma_m1P = gamma_func(s, l1, m1);
        gamma_mN = gamma_func(s, l, -m);
        gamma_mP = gamma_func(s, l, m);

        lambda_mN = lambda_func(s, l, -m);
        lambda_mP = lambda_func(s, l, m);

        u4 = -1i*p/2 * gamma_m1N .* sigma_mP .* sigma_l1N ...
            + 1i*p/2 * lambda_mP .* sigma_mP .* sigma_l ...
            + 1i*p/2 * gamma_mP .* sigma_mP .* sigma_l1P ...
            + 1i*p/2 * gamma_m1P .* sigma_mN .* sigma_l1N ...
            + 1i*p/2 * lambda_mN .* sigma_mN .* sigma_l ...
            - 1i*p/2 * gamma_mN .* sigma_mN .* sigma_l1P;

        u5 = -p/2 * gamma_m1N .* sigma_mP .* sigma_l1N ...
            + p/2 * lambda_mP .* sigma_mP .* sigma_l ...
            + p/2 * gamma_mP .* sigma_mP .* sigma_l1P ...
            - p/2 * gamma_m1P .* sigma_mN .* sigma_l1N ...
            - p/2 * lambda_mN .* sigma_mN .* sigma_l ...
            + p/2 * gamma_mN .* sigma_mN .* sigma_l1P;

        kappa1 = kappa_func(s, l1, m1);
        kappa = kappa_func(s, l, m);

        u6 = 1i*p * kappa1 .* sigma_m .* sigma_l1N ...
            + 1i*p * s*m/l/(l+1) .* sigma_m .* sigma_l ...
            + 1i*p * kappa .* sigma_m .* sigma_l1P;

        u = x(1) * u1 + x(2) * u2 + x(3) * u3 ...
            + x(4) * u4 + x(5) * u5 + x(6) * u6;
    end
end

function gamma = gamma_func(s, l, m)
%   gamma is defined as gamma_{l, m}^s
    if l ~= 0
        gamma = sqrt( ...
            (l^2-s^2) * (l-m) * (l-m-1) ...
                / l^2 / (2*l-1) / (2*l+1) ...
            );
    else
        gamma = 0;
    end
end

function lambda = lambda_func(s, l, m)
%   lambda is defined as lambda_{l, m}^s
    if l ~= 0
        lambda = s * sqrt((l-m) * (l+m+1)) ...
            / l / (l+1);
    else
        lambda = 0;
    end
end

function kappa = kappa_func(s, l, m)
%   kappa is defined as kappa_{l, m}^s
    
    if l ~= 0
        kappa = sqrt( ...
            (l^2-m^2) * (l^2-s^2) ...
                / l^2 / (2*l-1) / (2*l+1) ...
            );
    else
        kappa = 0;
    end
end