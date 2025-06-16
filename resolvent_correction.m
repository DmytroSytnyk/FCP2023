function [A, cz0] = resolvent_correction(fOp, z, z0, n, Varg)

%RESOLVENT_CORRECTION performs resolvent correction
%
%INPUT:
% fOp       - function that calculates the action of operator H on Varg.
%             It has prototype: @(x)
% z         - (scalar or vector) point to evaluate resolvent correction
%             sum_{j=1}^n (H-z0*I)^(j-1)/(z-z0)^j Varg
% z0        - correction point
% Varg_grid - space grid where the given argument is evaluated
% Varg      - vector value of the argument on the grid Varg_grid
% n         - order of resolvent correction

%OUTPUT:
% A         - the value of resolvent correction for Varg 
% cz0       - correction point
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
if (n==0)       % No correction is needed 
  A = zeros(size(Varg)); cz0 = z0;
else
  cz0 = z0;
  if isinf(z0)
    MOp_powers = op_power(fOp, n-1, Varg, true);
    cz0 = min(real(z)) - max(vecnorm(MOp_powers,Inf,1));
  elseif (z0 == 0)
    MOp_powers = op_power(fOp, n-1, Varg, true);
  end
  % If cz0 <> 0 then we need to calculate the powers of H-z0*I.
  if (cz0 ~= 0)
    fOpdiff = @(x) fOp(x) - cz0*x;
    MOp_powers = op_power(fOpdiff, n-1, Varg, true);
  end
  denom = @(x,p) 1./(x.^p);
  n_z = numel(z);
  zprod = bsxfun(denom, reshape(z-cz0,[1 n_z]), (1:n).');
  A = MOp_powers*zprod;
end
end

