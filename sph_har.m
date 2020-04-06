function Ylm = sph_har(phi, theta, l, m_index)
%SPH_HAR    spherical harmoics relative to IURs of SO(3)
%   Ylm = sph_har(phi, theta, l)
%   Ylm = sph_har(phi, theta, l, m_index)
%
%   Input:
%    - phi: azimuth, between 0 and 2*pi
%    - theta: elevation, between 0 and pi
%    - l: positive integer, dual index
%    - m: integer between -l and l, index of Ylm
%   Output:
%    - Ylm: spherical harmonics

%-- Auther: hshi17 11/17/18 --%

    if isvector(phi) && isvector(theta)
        [Phi, Theta] = meshgrid(phi, theta);
    else
        Phi = phi;
        Theta = theta;
    end
    
    if ~isequal(size(Phi), size(Theta))
        disp('size of phi, theta not equal, please input a grid of angles');
    end
    
    if l == 0
        Ylm = 1;
        return;
    end
    
    m = -l:1:l;
    
    % calculate P
    Plm = legendre(l,cos(theta)); % for m >= 0
    m_neg = -m(m<0).';
    Plm_neg = (-1).^m_neg .* ...
        factorial(l-m_neg)./factorial(l+m_neg) .* ...
        flip(Plm(2:end, :, :),1);
    Plm = cat(1, Plm_neg, Plm);
    
    % calculate Y
    a = (2*l+1)*factorial(l-m);
    b = 4*pi*factorial(l+m);
    c = sqrt(a./b);
    M = repmat(m.', 1, size(phi,1), size(phi,2));
    Phi = repmat(permute(phi, [3,1,2]), [length(m),1,1]);
    d = exp(1i*M.*Phi);
    Ylm = permute(c.' .* Plm .* d, [2,3,1]);
    
    if nargin == 4
        if abs(m_index) > l
            Ylm = 0;
        else
            Ylm = Ylm(:, :, m_index+l+1);
        end
    end

end
