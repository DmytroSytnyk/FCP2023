function [u,H,r] = resolvent_fd_1d_epdo_bc_dirichlet(z, grid, vv, eps, a, b, C,method,method_opts)

% RESOLVENT_FD_1D_EPDO_BC_DIRICHLET finited-difference solver of resolvent equation for elliptic partial differential operator in 1-D
%
% [u,H,r] = resolvent_fd_1d_epdo_bc_dirichlet(z,grid,vv,eps,a,b, C,method,method_opts)
% Numerically solve a resolvent equation
%       (z*I - H)u = vv
% where Hu = - eps * d2u/dx^2 + a * du/dx + b*u is the differential operator
% defined on the interval x \in (0, L) accompained by the following
% boundary conditions:
%       u(0,0) = C(1)
%       u(T,L) = C(2)
% using finite difference scheme on a user-defined grid.
% To make the solving procedure compatible with eigenvalue estimation
% technique we solve the equation
%       (H - zI)u = -vv
% instead of the former one.
%
%INPUT:
% z           - z for the resolvent
% grid        - finite-difference grid
% vv          - right-hand side evaluated on grid
% eps         - diffusion coefficient
% a           - velocity
% b           - source term
% C           - rhs coefficients for the BC (two element vector)
% method      - Matlab method to solve resolvent equation
% method_opts - options for the chosen method
%
%OUTPUT:
% u     - the solution of the BVP (a.k.a resolvent equation)
% H     - discretized version of the operator
% r     - residual for the solution of the discretized system
%
% For the finite-difference formulas see [1]. The handling of Neumann BC is
% adopted from [2].
%
%REFERENCES:
% 1.  Bowen Derivative formulas and errors for non-uniformly spaced points
% 2.  Long Cheng FINITE DIFFERENCE METHODS FOR POISSON EQUATION
% NOTE! There is an error in formula (A 3b) of [1] (factor 1/3 is missing 
% in front of each term) that is corrected here.
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin == 6)
  C = [0,0];
  method = '\';
  method_opts = {};
elseif (nargin == 7)
  method = '\';
  method_opts = {};
elseif (nargin == 8)
  method_opts = {};
elseif (nargin == 9)
  % Do nothing here
else
  error('Wrong number of input arguments');
end

if strcmpi(method,'\') || strcmpi(method,'mldivide')
  SOLVER = @mldivide;
elseif strcmpi(method,'gmres')
  SOLVER = @(A,b) gmres(A,b,method_opts{:});
elseif strcmpi(method,'lsqr')
  SOLVER = @(A,b) lsqr(A,b,method_opts{:});
else
  error('User-provided solution method is not supported');
end

if any(isnan(grid))
  error('Some of the gridpoints are NaN.');
elseif isempty(grid)
  error('Parameter grid is empty. Please provide a finite-diference grid.');
end

% if (C(1) ~= 0) || (C(2) ~= 0)
%   error('This method support Dirichlet BC only. For general BC consider using different resolvent implementation.');
% end

%% Function body
EPS = 1e-5;
n = numel(grid);
vgrid = reshape(grid,1, n);

if any(abs(C - [vv(1),vv(n)]) > EPS)
%   assert(1,fprintf("Provided argument and the boundary conditions are not compatible. \nThe numerical method won't converge in the sup-norm."));
end
if (numel(vv) > (1e3-1))
  Id = speye(n);
else
  Id = eye(n);
end
% right−hand side of Ay = vb
vb = zeros(n,1);

% Assemble the discretized version of -eps * d^2 u/dx^2 into A
dd = @(x,y,z) 2./((x-y).*(x-z));
x1 = vgrid(1:n-2);
x2 = vgrid(2:n-1);
x3 = vgrid(3:n);
dd1 = dd(x1,x2,x3);
dd2 = dd(x2,x1,x3);
dd3 = dd(x3,x1,x2);
if issparse(Id)
  Bin = [[dd1,0,0];[0,dd2,0];[0,0,dd3]].';
  H = spdiags(Bin,-1:1,n,n);
else
  H = diag([dd1 0],-1)+diag([0 dd3],1)+diag([0 dd2 0]);
end
H = -eps*H;

% Assemble the discretized version of a * du/dxd
if (a ~= 0)
  d1 = -0.5*dd1.*(x2+x3-2*x1);
  d2 = -0.5*dd2.*(x1+x3-2*x2);
  d3 = -0.5*dd3.*(x1+x2-2*x3);
  if issparse(Id)
    Bin = [[d1/3,0,0];[0,d2/3,0];[0,0,d3/3]].';
    H1 = spdiags(Bin,-1:1,n,n);
  else
    H1 = (diag([d1 0],-1)+diag([0 d3],1)+diag([0 d2 0]))/3;
  end
  H = H + a*H1;
end

% Account for b
H = H + b*Id;

%% Boundary values
% In this implementation of Dirichlet BC
% BC: u(t,0) = 0
% BC: u(t,L) = 0
% we factor out the boundary nodes. Hence the size of the matrix is n-2.

vb = vv(2:n-1);
vb = -vb;
H = H(2:n-1,2:n-1);

nz = numel(z);
u = zeros(n,nz);
if (nargout > 2)
  r = zeros(n,nz);
end
for j = 1:nz
  %   y = (H - z(j)*Id(2:n-1,2:n-1))\ vb;
  y = SOLVER(H - z(j)*Id(2:n-1,2:n-1),vb);
  %     y = lsqminnorm((A - z(j)*eye(n)),vb);
  u(2:n-1,j) = y;
  if (nargout > 2)
    r(2:n-1,j) = ((H - z(j)*Id(2:n-1,2:n-1))*y-vb);
  end
end

end
