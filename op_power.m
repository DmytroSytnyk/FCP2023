function A = op_power(fOp, n, Varg, all_powers)

%OP_POWER evaluates action of the power of operator on the given vector

%INPUT:
% fOp(Varg)  - function that calculates the action of operator H on Varg
% n          - power of operator that needs to be evaluated
% Varg       - the value of the argument on the grid Varg_grid
% all_powers - if true return the intermediate operator powers 0,1,..,n-1

%OUTPUT:
% A - result of the evaluation (might be an array if all_powers=1)
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
x = Varg;
if (all_powers)
  MOp_powers = zeros(numel(x), n+1);
  MOp_powers(:,1) = x;
end
for j=1:n
  x = fOp(x);
  if (all_powers)
    MOp_powers(:,j+1) = x;
  end
end
if (all_powers)
  A = MOp_powers;
  return;
end
A = x;
end
