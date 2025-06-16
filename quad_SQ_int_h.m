function [r] = quad_SQ_int_h(vf, vx, a, b, N, h)

%QUAD_SQ_INT_H approximates the integral of the function over [a,b] using Sinc Quadrature rule
%
%INPUT:
% vf - vector of function values (or matrix if the function is
%      column-vector valued)
% a  - left endpoint of the inteval
% b  - right endpoint of the interval
% N  - quadrature's discretization parameter (actual size of grid is 2*N+1)
% h  - the grid stepsize
%
%OUTPUT:
% r  - result of the quadrature (matrix with the size size(vf,1)xnumel(vx))

%% Parse input arguments
if (size(vf,2) ~= 2*N+1) || (size(vx,2) ~= 2*N +1)
  error('Input vector dimensions are incorrect');
end

%% Function body
% wi = wq
% wi = h*ones(size(vf))./map_axis2int_diffinverse(a,b,vx);
wq = h*map_axis2int_diffdirect(a,b,h*(-N:N));
r = sum(vf.*repmat(wq,size(vf,1),1),2);
end
