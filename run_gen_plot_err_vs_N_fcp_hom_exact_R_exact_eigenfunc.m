function [v_alpha,v_N, err] = run_gen_plot_err_vs_N_fcp_hom_exact_R_exact_eigenfunc()
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

environment

%% Error norms 
err_norm_x = Inf;
err_norm_t = Inf;

% Lists of parameters to perform the calculations on
% The following conde will be executed for each combination of the
% parameters from the list
% alpha = [0.1, 0.5, 0.7, 0.9, 1, 1.1, 1.5, 1.7, 1.9];
alpha = [1, 0.7, 1.2, 0.5, 1.5, 0.3, 1.7, 0.1, 1.9];
% alpha = [0.1, 0.4, 1, 1.4, 1.9];
%  alpha = [1, 0.7, 0.5, 0.3, 0.1];

N = [32:32:2048];
% N = power(2,3:9);

%% Grid parameters
Nt = 200;
Nx = 100;
% t = [0];
x = (0:1/Nx:1).';

%% Resolvent parameters
% We take first and second initial values to be eigenfunction of A that
% correpsonds to eigenvalues with indices ev_n0, ev_n1
ev_n0 = 1;              % Eigenfunction index for u_0
ev_n1 = 4;              % Eigenfunction index for u_1
op_c = 1e2;             % We consider A = -op_c*d^2u/dx^2

%% Script body
method = 'hyperbola';
alpha = unique(alpha,'stable');
% desc=sprintf("D= %.2f, u_0,u_1 - %d,%d-th eigenfunction.\n Contour: %s",op_c, ev_n0,ev_n1,method);
desc = [];
fig_opts = {'Position',  [10, 50, 520, 450]};
figure(fig_opts{:})
fig = gca();

T = 1; t = [0:1/Nt:T];

plot_opts = {'Marker','.', 'MarkerSize',30, 'LineStyle','-', 'LineWidth',3};
[~,fig] = plot_err_vs_n(@error_fcp_hom_exact_R_exact_eigenfunc,alpha(alpha<=1),N,desc,fig,plot_opts,t,x,err_norm_t,err_norm_x);

hold on

T = 5;t = [0:1/Nt:T];

plot_opts = {'Marker','x', 'MarkerSize',11, 'LineStyle','-', 'LineWidth',3};
[~,fig] = plot_err_vs_n(@error_fcp_hom_exact_R_exact_eigenfunc,alpha(alpha>1),N,desc,fig,plot_opts,t,x,err_norm_t,err_norm_x);
xlim([0,max(N)+min(N)])
ylim([1e-16,1])
h = legend(fig);
set(h,'fontsize',32)
set(gca,'FontSize',26)
xlabel(''); ylabel('');

function [err] = error_fcp_hom_exact_R_exact_eigenfunc(fcp_alpha,fcp_N,t,x)
% Flags (0,1) on whether to include the contribution from the first 
% (include_iv0 = 1) and the second (include_iv1 = 1) propagator
% This essentially allows to debug two propagators separately
include_iv0 = 1;
include_iv1 = 1;

%% Spectrum shape parameters
fcp_spectral_estimator = 'manual';
fcp_contour_vertex   = NaN;   % NaN - Auto
fcp_spectral_vertex  = 0;     % Overridden by zero for 'hyperbola' method
fcp_spectral_angle   = pi/60;
fcp_critical_vertex  = -pi/8; 
fcp_critical_angle   = pi/2;

%% Correction paramters
fcp_correction_order = 1;   %Only integer values 0,1 are implemented
fcp_correction_point = 0;   %NaN - Auto

fcp = struct('frac_order',fcp_alpha, 'N',fcp_N, ...
  'spectral_estimator',fcp_spectral_estimator, ...
  'contour_vertex',fcp_contour_vertex, ...
  'spectral_vertex',fcp_spectral_vertex, ...
  'spectral_angle',fcp_spectral_angle, ...
  'critical_vertex',fcp_critical_vertex, ...
  'critical_angle',fcp_critical_angle, ...
  'correction_order',fcp_correction_order, ...
  'correction_point',fcp_correction_point);

%% Resolvent of spatial operator and initial value for the problem 

% NOTE!!! 
% The resolvent function fRe should adhere to the prototype R = fRe(z,Varg): 
%             R = (z*I + A)^(-1)Varg
% where A is the operator with spectrum compatible with the spectrum 
% parameters fcp_spectral_vertex, fcp_spectral_angle defined above.

a = max(x);
u0_res = @(z,u) resolvent_exact_eigenfunc(z,ev_n0,x,a,op_c);
u1_res = @(z,u) resolvent_exact_eigenfunc(z,ev_n1,x,a,op_c);
[vu0,~] = eigenfunction(a,ev_n0,x);
[vu1,~] = eigenfunction(a,ev_n1,x);
if (include_iv0 == 0)
  vu0 = zeros(size(vu0));
  u0_res = @(z,u) zeros(numel(u),numel(z));
end
if (include_iv1 == 0)
  vu1 = zeros(size(vu1));
  u1_res = @(z,u) zeros(numel(u),numel(z));
end
ev0 = eigenvalue(a,ev_n0,op_c);
ev1 = eigenvalue(a,ev_n1,op_c);


u_appr = fcp_cauchy_hom_appr(t,fcp,u0_res,u1_res,vu0,vu1,method);
u_ex   = fcp_cauchy_hom_ex_sol_eigenfunction(t,fcp_alpha,ev0,vu0,ev1,vu1);
err    = u_ex-u_appr;
end

end
