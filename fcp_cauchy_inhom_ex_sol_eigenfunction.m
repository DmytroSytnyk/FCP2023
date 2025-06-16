function [fcp_prop_t] = fcp_cauchy_inhom_ex_sol_eigenfunction(t,alpha,grid,op_c,ev_f0,f0_c,ev_df,df_t,NI,NJ,d)

%FCP_CAUCHY_INHOM_EX_SOL_EIGENFUNCTION - exact solution of inhomogenoeous Fractional Cauchy problem with RHS having eigenfunction shape
%
% [fcp_prop_t] =
% frac_cauchy_inhom_ex_sol_eigenfunction(t,alpha,op_c,ev_f0,f0_c,ev_df,df_t,NI,NJ,df_strip_height)
% Evaluation solution of fractional Cauchy problem with initial values
% that are eigenfunctions of the spatial operator coefficient.
%
%INPUT:
% t     - evaluation time (row vector)
% alpha - order of fractional derivative
% ev_f0 - eigenvalue for the first intitial conditions
% f0_c  - the value of RHS at 0 f(0)
% ev_df - eigenvalue for the derivative of f part of the solution
% df_t  - the function df = @(t) that returns RHS derivative for a given t
% NI    - discretization parameter for the outer integral
% NJ    - discretization parameter for the Riemann-Liouville integral
% d     - analyticity strip half-height of df(map_axis2int_direct(t))
%
%OUTPUT:
% fcp_prop_t - array of solution vector at each timestep of the size 
%              numel(grid) x numel(t)
%REFERENCES:
% [1] Gorenflo, R. et al. 2020. Mittag-Leffler Functions, Related Topics
% and Applications. Springer.
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input arguments
% Check whether df_t depents on t or not 
if ( size(df_t(t(1))) == size(df_t(t)) )
  df_t = @(t) repmat(df_t(t(1)), size(t));
end

%% Function body
t = reshape(t,1,[]);
grid = reshape(grid, [],1);
a = max(grid);

% First term: J_alpha * S * f(0)
% The variable T is used within RLI_approx only to check the sanity of t 
T = max(t)+1;
f0_eig = eigenvalue(a,ev_f0,op_c);
df_eig = eigenvalue(a,ev_df,op_c);
% S_f0 = @(t) ml(-f0_eig*t.^alpha,alpha,1).*f0_c;
% fcp_prop_t = RLI_approx(S_f0, alpha, t, 0, T, NJ, pi/4).*eigenfunction (a, ev_f0, grid);
% The above RL integral of S_f0 can be evaluated explicitely [1, Eq. 3.7.44]
fcp_prop_t = (1 - ml(-f0_eig*t.^alpha,alpha,1)).*(f0_c/f0_eig).*eigenfunction (a, ev_f0, grid);
% Second term: int(S * J_alpha * f'(s), s = 0..t)
% Because of the substitution s = t*s, made prior to the application of 
% the sinc-quadrature, the latter one is essentially always applied on [0,1]
fcp_prop_t = fcp_prop_t + ((t(:).*quad_SQ_int(@S_J_df,0,1,NI,d)).').*eigenfunction (a, ev_df, grid);

%% Local functions
function r = S_J_df(s)
  S_J_int = @(t,s) ml(-df_eig*(t(:).*(1-s)).^alpha, alpha, 1).* ...
            RLI_approx(df_t, alpha, t(:).*s, 0, T, NJ, d);
  r = bsxfun(S_J_int,t(:),s);
end
% function r = S_J_int(t,s) 
%   r = t.*ml(-op_c*t.*(1-s).^alpha,alpha,1).* ...
%   RLI_approx(df_t, alpha, s, 0, 1, NJ, df_strip_height);
% end

end
