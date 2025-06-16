function [u,k_out,grid_out,a_out,op_c_out] = resolvent_exact_eigenfunc(z, k, grid, a, op_c, arg)

%RESOLVENT_EXACT_EIGENFUNC Evaluate 1-D Laplacian resolvent with Dirichlet BC on its eigenvector
%
% [u,k_out,grid_out,a_out,op_c_out] = resolvent_exact_eigenfunc (z, k, grid, a, op_c,s) 
% Evaluate the action of (z*I - A)^(-1), where
%             A = -op_c*d^2/dx^2 
% is defined on the interval x \in [0,a] with zero BC,
% on its eigenfunction
%             v_k = sin(pi/a*k*x)
% the corresponding eigenvalue of A is  op_c*(pi/a*k)^2
%
%INPUT:
% z    - argument of the resolvent
% k    - eigenfunction number
% grid - grid on the interval [0, a]
% a    - length of the interval 
% op_c - diffusion coefficient
% arg  - argument (or the scalar value that is multiplied by the result)
%
%OUTPUT:
% u    - resulting function discretized on the grid 
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin == 2)
  a = 1;
  grid = (0:1/100:a).';
  op_c = 1;
  arg = 1;
elseif (nargin == 3)
  a = max(grid(:));
  op_c = 1;
  arg = 1;
elseif (nargin == 4)
  op_c = 1;
  arg = 1;
elseif (nargin == 5)
  arg = 1;
elseif (nargin == 6)
  % Do nothing
else
  error('Wrong number of input arguments.');
end
grid = reshape(grid,[],1);

%% Function body
ev = eigenvalue(a,k,op_c);
ef = eigenfunction(a,k,grid);
if isequal(size(arg), size(grid))
  u = arg./(z - ev);
elseif (numel(arg) == 1)
  u = arg.*ef./(z - ev);
else
  error('Size of the provided argument s is incompatible with grid.');
end
% Optionally output resolvent parameters when needed
if (nargout > 1)
  k_out = k;
  grid_out = grid;
  a_out = a;
  op_c_out = op_c;
end

end

function [f] = efunc(a,k,x)
  f = eigenfunction (a, k, x);
end
