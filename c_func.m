function c = c_func (l, n)
%C_FUNC     a simple function used in diff_iur_so3.m and diff_iur_se3.m
%           c is defined as c_n^l

%-- Auther: hshi17 12/20/18 --%

    c = sqrt((l-n) .* (l+n+1)) .* (n<=l);
end