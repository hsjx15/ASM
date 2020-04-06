function u = diff_iur_so3(X, lambda, m, n)
%DIFF_IUR_SO3   Differentiation of IUR for SO(3)
%   u = DIFF_IUR_SO3(X, lambda)
%   u = DIFF_IUR_SO3(X, lambda, m, n)
%
%   Input:
%    - X: 3*3 matrix in so(3)
%    - lambda: positive integer, dual index
%    - m, n: integers between -lambda and lambda, indices of u
%   Output:
%    - u: differentiation of IUR for SO(3), lambda*lambda or m*n matrix

%-- Auther: hshi17 11/17/18 --%

    if nargin == 2
        m = -lambda:1:lambda;
        n = m;
    end

    [N, M] = meshgrid(n, m);   
    x = so32vec(X);

    c_n = c_func(lambda, -N);
    cn = c_func(lambda, N);
    
    u1 = -1i/2 * c_n .* ((M+1)==N) ...
        -1i/2 * cn .* ((M-1)==N);
    u2 = 1/2 * c_n .* ((M+1)==N) ...
        -1/2 * cn .*((M-1)==N);
    u3 = -1i * N .* (M==N);

    u = x(1) * u1 + x(2) * u2 + x(3) * u3;  
end