function [cont, cont_deriv] = cont_hyperbolic (x, a0, a_i, b_i)

% CONT_HYPERBOLIC parametrization of hyperbolic contour
% x \in (-\infty; \infty)
%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved. 

cont = a0 + a_i*cosh(x) - 1i * b_i*sinh(x);
cont_deriv = a_i*sinh(x) - 1i * b_i*cosh(x);
end

