function U = IUR_SO3(R, lambda, m, n)
%IUR_SO3    Irreducible Unitary Representations for SO(3)
%   U = IUR_SO3(R, lambda)
%   U = IUR_SO3(R, lambda, m, n)
%
%   Input:
%    - R: 3*3 matrix in SO(3)
%    - lambda: positive integer, dual index
%    - m, n: integers between -lambda and lambda, indices of U
%   Output:
%    - U: IUR for SO(3), lambda*lambda or m*n matrix

%-- Auther: hshi17 11/17/18 --%
    
    if nargin == 2
        m = -lambda:1:lambda;
        n = m;
    end
    
    dphi = pi/1e4;  % integral segment
    phi = 0:dphi:2*pi;
    [N, M, Phi] = meshgrid(n, m, phi);

    % derive alpha, beta, gamma angles from G:
    [alpha, beta, gamma] = SO3param(R);

    % get the integral part of the wigner-d function:
    cb = cos(beta/2);
    sb = sin(beta/2);
    ip = 1i*Phi/2;
    wigner_d_int = sum( ...
        (cb * exp(ip) + 1i * sb * exp(-ip)) .^ (lambda-N) ...
        .* (cb * exp(-ip) + 1i * sb * exp(ip)) .^ (lambda+N) ...
        .* exp(1i * M .* Phi) * dphi, 3);

    % calculate wigner-d function:
    [N, M] = meshgrid(n, m);
    wigner_d =  (1i) .^ (M-N) / 2 / pi ...
        .* sqrt(factorial(lambda-M) .* factorial(lambda+M) ...
        ./ factorial(lambda-N) ./ factorial(lambda+N)) ...
        .* wigner_d_int;

    % calculate U:
    U = exp(-1i*M*alpha) .* wigner_d .* exp(-1i*N*gamma);
end