function h = quad_SQ_get_h(N, d, alpha)

%QUAD_SQ_GET_H calculates quadrature step-size
% N     - discretization parameter for sinc-quadrature formula 
% d     - height of the analyticity strip of needed for the Sinc-quadrature
% alpha - exponential decay order of the integrand |f(t)| < C e^(-alpha|t|)
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input arguments
if (nargin == 2)
  alpha = 1;
end

%% Function body
h = sqrt(2*pi*d/(N*alpha));
end

