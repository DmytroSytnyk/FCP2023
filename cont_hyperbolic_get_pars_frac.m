function [a0,ai,bi,d] = cont_hyperbolic_get_pars_frac(alpha,sa0,sphi,ca0,cphi,d_in)

% CONT_HYPERBOLIC_GET_PARS calculate parameters of the sector-enveloping hyperbola 
%
% [a0,ai,bi,d] = cont_hyperbolic_get_pars_frac(alpha,sa0,sphi,ca0,cphi) 
% Calculate the parameters of hyperbolic contour
%     z = a0 + a_i*cosh(x) - 1i * b_i*sinh(x),  x \in (-\infty; \infty)
% and the width of the strip d based on two sectors: spectral and critical.
%
%INPUT:
% alpha - fractional order parameter
% sa0   - abscissa of the spectral sector's vertex 
% sphi  - angle of the spectral sector 
% ca0   - abscissa of the critical sector's vertex
% cphi  - angle of the critical sector 
% d_in  - manual override for the values of d
%
%OUTPUT:
% a0 - contour parameter (see above)
% ai - positive contour parameter (see above)
% bi - positive contour parameter (see above)
% d - height of the strip of analyticity needed for the Sinc-quadrature
%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved. 


%% Parse input parameters
if (nargin == 3)
  ca0 = -pi/2;
  cphi = pi/2;
  d_in = pi/2;
elseif (nargin == 4)
  cphi = pi/2;
  d_in = pi/2;
elseif (nargin == 5)
  d_in = pi/2;
elseif (nargin == 6)
  % Do nothing
else
  error('Wrong number of input parameters');
end

if (sa0 < 0)
  error('The case sa0 < 0 is not implemented. Set sa0 = 0.');
end
if (alpha<=0)||(alpha>=2)
  error('The fractional parameter alpha should be within (0, 2).');
end


%% Function body
sphia = min(pi,(pi-sphi)/alpha);
sa0a = sa0^(1/alpha);
d = (sphia+cphi-pi)/2;
if (d_in < d) && (d_in > 0)
  d = d_in;
  warning('The value of d is manually overriden! Make sure you know what you are doing.')
end
a0 = ca0;
ai = (sa0a-a0)*(cos(d) + sin(d)*tan(sphia));
bi = (sa0a-a0)*(sin(d) - cos(d)*tan(sphia));
if (ai <= 0)|| (bi <= 0)
  error('Internal error. The contour fails to encircle the spectrum');
end
end
