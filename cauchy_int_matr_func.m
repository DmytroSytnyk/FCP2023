function [VOpFunc, res_array, new_cor_point] = cauchy_int_matr_func (fOp, fRes, fCont, Mfunc, Varg, cor_order, cor_point, quad_h, res_array_in, quad_idx)

%CAUCHY_INT_MATR_FUNC approximate Cauchy integral.
%
%   cauchy_int_matr_func(fOp, fRes, fCont, Mfunc, Varg_grid, Varg, 
%   cor_order, cor_point, quad_h, quad_d, res_array_in, quad_idx)  
%   Calculates approximation of the Cauchy integral using trapezoidal
%   rule with a given step-size quad_h and number of grid points equal to
%   the number of columns of Mfunc.
%   
%INPUT:
%
% fOp          - function @(Varg) that calculates the action of operator 
% fRes         - function @(z, Varg) that calculates action of resolvent 
% fCont        - function @(x) that return the value on a contour and 
%                its derivative at x
% Mfunc        - matrix representation of the operator function 
%                where f(t_k,fCont(quad_h * k)) for a list of values t_k
%                and k goes through quad_idx 
% Varg         - the value of the argument (column vector) 
% cor_order    - order of resolvent correction
% cor_point    - point for the correction
% quad_h       - step size of contour quadrature 
% quad_d       - strip height parameter of contour quadrature
% res_array_in - pre-calculated vector of corrected resolvents (optional)
% quad_idx     - vector indices to evaluate the quadrature on 
%                (by default equal to (-N:N))
%
%OUTPUT:
%
% VOpFunc      - result of approximation
% res_array    - matrix of calculated resolvent values for a given quad_idx 
%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
PARALLELIZE_RESOLVENT_EVALUATIONS=false;

N = floor(size(Mfunc,2)/2);
if isempty(quad_idx) 
  idx = (-N:N); 
else 
  idx = quad_idx; 
end
n_idx = numel(idx);
n_arg = numel(Varg);
[cont_quad_points, cont_deriv] = fCont(quad_h * idx);
if ~(all(isfinite(cont_quad_points)) && all(isfinite(cont_deriv)))
  error('Some of the contour or derivative quadrature points are not finite');
end
if isequal(size(res_array_in),[n_arg, n_idx])
  % Reuse the input array of resolvents
  res_array = res_array_in;
else
  % Calculate resolvents and corrections
  if PARALLELIZE_RESOLVENT_EVALUATIONS && (~isempty(gcp('nocreate')))
    MRes = zeros(numel(Varg),numel(cont_quad_points));
    parfor i=1:numel(cont_quad_points)
      MRes(:,i) = fRes(cont_quad_points(i), Varg);
    end
  else
    MRes = fRes(cont_quad_points, Varg);
  end
  [MRes_cor,new_cor_point] = resolvent_correction(fOp, cont_quad_points, cor_point, cor_order, Varg);
  res_array = MRes - MRes_cor;
end
VOpFunc = quad_h/2/pi/1i*(res_array * bsxfun(@times, Mfunc, cont_deriv).'); 
end

