function [r] = quad_SQ_int(f, a, b, N, d)

%QUAD_SQ_INT approximates the integral of the function over [a,b] using Sinc Quadrature rule
%
% [r] = quad_SQ_int(f, a, b, N, d)
%
%INPUT:
% f - a vectorizable function y = @(x) that return values of integrand
% a - left endpoint of the inteval 
% b - right endpoint of the interval
% N - discretization parameter for sinc-quadrature formula 
% d - height of the analyticity strip of needed for the Sinc-quadrature
%
%OUTPUT:
% r  - result of the quadrature (matrix with the size size(vf,1)xnumel(vx))
%
%REFERENCES:
% [1] Stenger, F., 2012. Numerical methods based on Sinc and analytic 
% functions (Vol. 20). Springer Science & Business Media.
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input arguments
if (d<=0) || (N<0)
  error('Input parameters are out of range (d>0, N > 0)');
end

%% Function body
h = quad_SQ_get_h(N,d);
x = h*(-N:N);
X = map_axis2int_direct(a,b,x);
% Handle the case when f is a constant
if ( size(f(X(1))) == size(f(X)) ) 
  vf = repmat(f(X(1)),size(X));
else
  vf = f(X);
end
if (a == b)
  r = zeros(size(vf,1),1);
else
  r = quad_SQ_int_h(vf,X,a,b,N,h);
end
end
