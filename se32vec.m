function x = se32vec(X)
%SE32VEC
%
%   Input:
%   - X: 4*4 matrix in se(3)
%   Output:
%   - x: 6*1 vector

%-- Auther: hshi17 12/20/18 --%

    x = [X(3,2); X(1,3); X(2,1); X(1,4); X(2,4); X(3,4)];

end