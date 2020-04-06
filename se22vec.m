function x = se22vec(X)
%SE22VEC
%
%   Input:
%    - X: 3*3 matrix in se(2)
%   Output:
%    - x: 3*1 vector

%-- Auther: hshi17 11/17/18 --%

    x = [X(1,3); X(2,3); X(2,1)];

end