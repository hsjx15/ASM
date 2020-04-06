function U = IUR_SE3(a, A, p, s, l1, l, m1, m)
%IUR_SE3    Irreducible Unitary Representations for SE(3)
%   U = IUR_SE3(a, A, p, s, l1, l)
%   U = IUR_SE3(a, A, p, s, l1, l, m1, m)
%
%   Input:
%    - a: translation, a 3*1 matrix
%    - A: rotation, a 3*3 matrix in SO(3)
%    - p: positive real number, first dual index
%    - s: integer, second dual index
%    - l1, m1, l, m: integers, sequenced indices of U, where
%           l1, l \ge abs(s), 
%           m1 \in {-l1, l1}, 
%           m \in {-l, l}
%   Output:
%    - U: IUR for SE(3), 11*m1*l2*m2 matrix

%-- Auther: hshi17 11/17/18 --%
    
    if ~isequal(size(a), [3,1])
        disp('a has wrong dimentsion, set a = [0;0;0]');
        a = zeros(3,1);
    end
    
    if ~isequal(size(A), [3,3])
        disp('A has wrong dimentsion, set A = eye(3)');
        A = eye(3);
    end
    
    if l1(1) < abs(s)
        disp('value of l1 is incorrect');
        U = NaN; return;
    end
    if l(1) < abs(s)
        disp('value of l is incorrect');
        U = NaN; return;
    end
    
    if nargin ~= 6 && nargin ~= 8
        disp('input error');
        U = NaN; return;
    end

    L = max(l1(end), l(end));
    U_size = (L+1)^2 - s^2;     % size for 2-dim matrix
    U = zeros(U_size, U_size);

    for i = 1:U_size
        for j = 1:U_size
            lc = ceil(sqrt([i j]+s^2))-1;   % current l1 and l
            
            if (lc(1) <= l1(end)) && (lc(2) <= l(end))    % in boundary
                mc = [i j] - lc .* (lc+1) -1 + s^2;     % current m1 and m

                U_temp = zeros(2*lc(2)+1,1);
                for k = -lc(2):1:lc(2)
                    U_temp(k+lc(2)+1) = ...
    IUR_SE3_trans(a, lc(1), mc(1), p, s, lc(2), k) ...
    * IUR_SO3(A, lc(2),k, mc(2)).';
                end
                U(i,j) = sum(U_temp);
            end
        end
    end
    
    if nargin == 8  % l1, l, m1, m are index numbers
        if m1 > l1
            disp('value of m1 is incorrect');
            U = 0; return;
        end
        if m > l
            disp('value of m is incorrect');
            U = 0; return;
        end
        U = U(l1*(l1+1)+m1-s^2+1, l*(l+1)+m-s^2+1);
    end
    
end

function U = IUR_SE3_trans(a, l1, m1, p, s, l, m)

    [az, el, r] = cart2sph(a(1), a(2), a(3));
    
    U = zeros((l1+l)-abs(l1-l)+1,1);
    
    for k = abs(l1-l):1:(l1+1)
        U(k-abs(l1-l)+1) = ...
            (1i)^k * sqrt((2*l1+1)*(2*k+1)/(2*l+1)) ...
            * besselj_sph(k, p*r) ...
            * ClebshGordan(k, 0, l1, s, l, s) ...
            * ClebshGordan(k, m-m1, l1, m1, l, m) ...
            * sph_har(az, el, k, m-m1);
    end
    
    U = sqrt(4*pi) * sum(U);
end

function J = besselj_sph(nu, z)
   J = sqrt(pi/2./z) .* besselj(nu+0.5, z);
end

function C = ClebshGordan(j1, m1, j2, m2, j, m)
    C = (-1).^(m+j1-j2) .* sqrt(2*j+1) ...
        .* Wigner3j([j1, j2, j], [m1, m2, -m]);
end