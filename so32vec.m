function x = so32vec(X)
%SO32VEC
%
%   Input:
%    - X: 3*3 matrix in so(3)
%   Output:
%    - x: 3*1 vector

%-- Auther: hshi17 11/17/18 --%

    x = [X(3,2); X(1,3); X(2,1)];

end