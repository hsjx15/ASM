function U = IUR_SE2(g, p, m, n)
%IUR_SE2    Irreducible Unitary Representions for SE(2)
%   IUR for SE(2) is an infinite dimensional matrix, the result here is a 
%   m x n dimensional approximation.
%
%   U = IUR_SE2(g, p, m, n)
%
%   Input:
%    - g: 3*3 matrix in SE(2)
%    - p: positive real number, dual index
%    - m, n: integers, indices of U
%   Output:
%    - U: IUR for SE(2), m*n matrix

%-- Auther: hshi17 11/17/18 --%

    [N, M] = meshgrid(n, m);
    
    [r, phi, theta] = SE2param(g);
    
    U = (1i).^(N-M) ...
        .* exp(-1i*(N*theta+(M-N)*phi)) ...
        .* besselj(N-M, p*r);   

end