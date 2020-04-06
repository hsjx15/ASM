function u = diff_iur_se2(X, p, m, n)
%DIFF_IUR_SE2    Differentiation of IUR for SE(2)
%   IUR for SE(2) is an infinite dimensional matrix, the result here is a 
%   m x n dimensional approximation.
%
%   u = DIFF_IUR_SE2(g, p, m, n)
%
%   Input:
%    - X: 3*3 matrix in se(2)
%    - p: positive real number, dual index
%    - m, n: integers, indices of u
%   Output:
%    - u: differentiation of IUR for SE(2), m*n matrix

%-- Auther: hshi17 11/17/18 --%

    [N, M] = meshgrid(n, m);
    
    x = se22vec(X);
    
    u1 = 1i*p/2 * (double(M==(N+1)) + double(M==(N-1)));
    u2 = p/2 * (double(M==(N+1)) - double(M==(N-1)));
    u3 = -1i * M .* double(M==N);
    
    u = x(1) * u1 + x(2) * u2 + x(3) * u3;

end