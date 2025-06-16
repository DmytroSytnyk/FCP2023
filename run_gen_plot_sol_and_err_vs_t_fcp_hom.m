function [] = run_gen_plot_sol_and_err_vs_t_fcp_hom()
%
% Copyright (C) 2022 - 2025, Dmytro Sytnyk. All rights reserved. 

environment
% Lists of parameters to perform the calculations on
% The following code will be executed for each combination of the
% parameters from the list
% alpha = [0.1, 0.5, 0.7, 0.9, 1, 1.1, 1.5, 1.7, 1.9];
alpha = [1, 0.7, 0.5, 0.3, 0.1];
% alpha = [1, 1.2, 1.5, 1.7, 1.9];

% N = [32];
% N = [32,64,128,256];
N = [8,16,32,64,128,256,512,1024];

%% Grid parameters
T = 1;
Nt = 100;
Nx = 100;
t = [0:T/Nt:T];
% t = [0];
x = (0:1/Nx:1).';
x_plot = 0.5;
%% Resolvent parameters
% We take first and second initial values to be eigenfunction of A that
% correpsonds to eigenvalues with indices ev_n0, ev_n1
ev_n0 = 1;              % Eigenfunction index for u_0
ev_n1 = 4;              % Eigenfunction index for u_1
op_c = 1;               % Scaler here we consider A = -op_c*d^2u/dx^2
% Flags (0,1) on whether to include the contribution from the first
% (include_iv0 = 1) and the second (include_iv1 = 1) propagator
% This essentially allows to debug two propagators separately
include_iv0 = 1;
include_iv1 = 1;

%% Plot exact solutions
desc =[sprintf("Exact solution $$a= %1.2f$$:", op_c); sprintf("$$u_0 = \\sin{(%d \\pi x)}$$, $$u_1 = \\sin{(%d \\pi x)}$$",ev_n0,ev_n1)];

fig_opts = {'Position',  [10, 50, 520, 450]};
figure(fig_opts{:})
fig = gca();
plot_opts = {'Marker','none', 'MarkerSize',30, 'LineStyle','-', 'LineWidth',4};
if include_iv0
  for i=find(alpha<=1)
    l_desc = sprintf('$\\alpha = %0.2f$',alpha(i));
    u_ex = @(t,x) sol_ex_fcp_hom_R_exact_eigenfunc(alpha(i),t,x);
    [~,fig] = plot_sol_vs_t(u_ex,x_plot,l_desc,fig,plot_opts,t,x);
  end
end
for i=find(alpha>1)
  l_desc = sprintf('$\\alpha = %0.2f$',alpha(i));
  u_ex = @(t,x) sol_ex_fcp_hom_R_exact_eigenfunc(alpha(i),t,x);
  [~,fig] = plot_sol_vs_t(u_ex,x_plot,l_desc,fig,plot_opts,t,x);
end
h = legend(fig);
xlabel('');
set(h,'fontsize',42)
set(gca,'FontSize',39)
title(desc,'interpreter','latex')

%% Plot approximations
l_desc = @(alpha)sprintf('$\\alpha = %0.2f$',alpha);
for j=1:numel(N)
  desc =[];
  desc =[sprintf("Approximate solution with $$N=%d$$:", N(j)); sprintf("$$u_0 = \\sin{(%d \\pi x)}$$, $$u_1 = \\sin{(%d \\pi x)}$$",ev_n0,ev_n1)];
  fig_opts = {'Position',  [10, 50, 520, 450]};
  figure(fig_opts{:})
  fig = gca();
  if include_iv0
    for i=find(alpha<=1)
      u_appr = @(t,x) sol_appr_fcp_hom_R_exact_eigenfunc(alpha(i),N(j),t,x);
      [~,fig] = plot_sol_vs_t(u_appr,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    end
  end
  for i=find(alpha>1)
    u_appr = @(t,x) sol_appr_fcp_hom_R_exact_eigenfunc(alpha(i),N(j),t,x);
    [~,fig] = plot_sol_vs_t(u_appr,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
  end
  h = legend(fig);
  xlabel('')
  set(h,'fontsize',42)
  set(gca,'FontSize',39)
  title(desc,'interpreter','latex')
end

%% Plot errors of approximation
l_desc = @(alpha)sprintf('$\\alpha = %0.2f$',alpha);
desc =[];
for j=1:numel(N)
  desc =[sprintf("Error of approximate solution with $$N=%d$$:", N(j)); sprintf("$$u_0 = \\sin{(%d \\pi x)}$$, $$u_1 = \\sin{(%d \\pi x)}$$",ev_n0,ev_n1)];
  fig_opts = {'Position',  [10, 50, 520, 450]};
  figure(fig_opts{:})
  fig = gca();
  if include_iv0
    for i=find(alpha<=1)
      u_err = @(t,x) abs(sol_ex_fcp_hom_R_exact_eigenfunc(alpha(i),t,x) - sol_appr_fcp_hom_R_exact_eigenfunc(alpha(i),N(j),t,x));
      [~,fig] = plot_sol_vs_t(u_err,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
      title(desc)
    end
  end
  for i=find(alpha>1)
    u_err = @(t,x) abs(sol_ex_fcp_hom_R_exact_eigenfunc(alpha(i),t,x) - sol_appr_fcp_hom_R_exact_eigenfunc(alpha(i),N(j),t,x));
    [~,fig] = plot_sol_vs_t(u_err,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
  end
  set(gca, 'YScale', 'log')
  h = legend(fig);
  xlabel('')
  set(h,'fontsize',42)
  set(gca,'FontSize',39)
  % ylim([1e-6,1])
  title(desc,'interpreter','latex')
end

  function [u_sol] = sol_ex_fcp_hom_R_exact_eigenfunc(fcp_alpha,t,x)
    a = max(x);
    ev0 = eigenvalue(a,ev_n0,op_c);
    ev1 = eigenvalue(a,ev_n1,op_c);
    [vu0,~] = eigenfunction(a,ev_n0,x);
    [vu1,~] = eigenfunction(a,ev_n1,x);
    if ~include_iv0
      % We need to make sure that both initial value and its resolven is zero
      vu0 = zeros(size(vu0));
    end
    if ~include_iv1
      vu1 = zeros(size(vu1));
    end
    u_sol   = fcp_cauchy_hom_ex_sol_eigenfunction(t,fcp_alpha,ev0,vu0,ev1,vu1);
  end
  function [u_appr] = sol_appr_fcp_hom_R_exact_eigenfunc(fcp_alpha,fcp_N,t,x)
    method = 'hyperbola';

    %% Spectrum shape parameters
    fcp_spectral_estimator = 'manual';

    fcp_contour_vertex   = NaN; % NaN - Auto
    fcp_spectral_vertex  = 0;   % Overridden by zero for 'hyperbola' method
    fcp_spectral_angle  = pi/60;
    fcp_critical_vertex = -pi/8;
    fcp_critical_angle  = pi/2;

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

    %% Resolvent of the spatial operator and initial value for the problem
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
    u_appr = fcp_cauchy_hom_appr(t,fcp,u0_res,u1_res,vu0,vu1,method);
  end
end