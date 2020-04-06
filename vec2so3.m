function X = vec2so3(x)
%VEC2SO3
%
%   Input:
%    - x: 3*1  or 1*3 vector
%   Output:
%    - X: 3*3 matrix in so(3)

%-- Auther: hshi17 11/17/18 --%

    X = NaN;

    if ~(isequal(size(x), [3,1]) || isequal(size(x), [1, 3]))
        disp('wrong input');
        return;
    end
    
    X = [0, -x(3), x(2);
         x(3), 0, -x(1);
         -x(2), x(1), 0];

end