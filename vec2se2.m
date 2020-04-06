function X = vec2se2(x)
%VEC2SE2
%
%   Input:
%    - x: 3*1  or 1*3 vector
%   Output:
%    - X: 3*3 matrix in se(2)

%-- Auther: hshi17 11/17/18 --%

    X = NaN;

    if ~(isequal(size(x), [3,1]) || isequal(size(x), [1, 3]))
        disp('wrong input');
        return;
    end
    
    X = [0, -x(3), x(1);
         x(3), 0, x(2);
         0, 0, 0];

end