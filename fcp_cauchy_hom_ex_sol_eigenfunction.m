function [fcp_prop_t] = fcp_cauchy_hom_ex_sol_eigenfunction(t,alpha,ev0,ef0,ev1,ef1)

% Evaluate solution of fractional Cauchy problem with initial values
% that are eigenfunctions of the spatial operator coefficient.
%
%INPUT:
% t     - evaluation time (row vector)
% alpha - order of fractional derivative
% ev0   - eigenvalue for the first initial conditions
% ef0   - eigenfunction that corresponds to ev0 evaluated at some grid 
%         (column vector)
% ev1   - eigenvalue for the first initial conditions
% ef1   - eigenfunction that corresponds to ev0 evaluated at some grid 
%         (column vector)
%
%OUTPUT:
% fcp_prop_t - matrix of discretized solution representations at given times
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
t = reshape(t,1,[]);

fcp_prop_t = ml(-ev0*t.^alpha,alpha,1).*ef0(:);
if (alpha > 1)
  fcp_prop_t = fcp_prop_t + ml(-ev1*t.^alpha,alpha,2).*(t.*ef1(:));
end
end
