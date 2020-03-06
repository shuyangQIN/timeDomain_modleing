function source = smooth_source(source,window_type)
%   - window type. Supported values are
%
%   'Bartlett'
%   'Bartlett-Hanning'   
%   'Blackman'
%   'Blackman-Harris'
%   'Blackman-Nuttall'
%   'Cosine'
%   'Flattop'
%   'Gaussian'
%   'HalfBand'
%   'Hamming'
%   'Hanning'
%   'Kaiser'
%   'Lanczos'
%   'Nuttall'
%   'Rectangular'
%   'Triangular'
%   'Tukey'
%

dim = size(source);
win = getWin(dim, window_type, 'Rotation', true , 'Symmetric', true, 'Param', 0.1);
source = real(ifftn(fftn(source) .* ifftshift(win)));

end