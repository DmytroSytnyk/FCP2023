function [] = run_gen_plot_sol_and_err_vs_t_fcp_cauchy__R_FD_Dirichlet()
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

environment
%% matlab cleanup and format
% clear all
% close all


% Lists of parameters to perform the calculations on
% The following code will be executed for each combination of the
% parameters from the list
% alpha = [0.1, 0.5, 0.7, 0.9, 1, 1.1, 1.5, 1.7, 1.9];
% alpha = [0.1, 0.4, 1, 1.4, 1.9];
% alpha = [0.5, 1];
% alpha = [1, 0.7, 0.5, 0.3, 0.1];
alpha = [1, 1.2, 1.5, 1.7, 1.9];

%  N = [256];
N = [16,32,64,128,256];
% N = [8,16,32,64,128,256,512,1024];

%% Grid parameters
T = 1;
Nt = 200;
Nx = 1e2;
t = [0:T/Nt:T];
% t = [0];
x = (0:1/Nx:1).';
x_plot = 0.1;

%% Error norm
err_norm_x = Inf;

%% Exact solution parameters
u_exact_f = @(t,x) x.^2 .* (x - t .^ 2 + 0.1e1 / 0.2e1) .* (x - 0.1e1);

%% Problem parameters
% Operator parameters
op_c = 1;

% Initial values
u0 = u_exact_f(0,x);
u1 = zeros(size(x));

% Right-hand side
% The corresponding functions f(0,x) and f'(t,x) are initialized below
f0_f = @(x) -12 * x .^ 2 + 3 * x + 1;
f0 = f0_f(x);
df_strip_height = pi/4; % The eye-shaped domain analyticity angle for df

%% Resolvent of spatial operator and the right hand side of the problem

% NOTE!!!
% The resolvent function fRe should adhere to the prototype R = fRe(z,Varg):
%             R = (z*I + A)^(-1)Varg
% where A is the operator with spectrum compatible with the spectrum
% parameters fcp_spectral_vertex, fcp_spectral_angle defined above.

% We take first and second initial values to be eigenfunction of A that
% (z, grid, vv, eps, a, b, C,method,method_opts)
a = max(x);
u0_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);
u1_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);
f0_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);
df_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);


fprintf('|A u_0| = %1.3e , |A df| = %1.3e .\n',vecnorm(u0_res(0,u0),err_norm_x,1),vecnorm(df_res(0,df_f(1,x,0.5)),err_norm_x,1));

%% Output file name
file_prefix = 'run_gen_plot_sol_and_err_vs_t_fcp_cauchy_R_FD_Dirichlet_x2_';
salpha = ['alpha_',sprintf('%1.1f_', sort(alpha))];
sN = ['N_',sprintf('%d_', sort(N))];
st = sprintf('t_%1.0f_%1.0f_',min(t),max(t));
sev = sprintf('Nt_%d_Nx_%d',Nt,Nx);
filename = [file_prefix,salpha,sN,st,sev,'.mat'];

%% Load data file
Lu_ex = NaN(numel(alpha),numel(x),numel(t));
Lu_appr = NaN(numel(N),numel(alpha),numel(x),numel(t));
storage =[];
if isfile(filename)
  fprintf("Associated data file:\n %s \n found.\n Loading ... ",filename);
  load(filename)
  fprintf("success.\n");
end

%% Plot exact solutions
desc =[sprintf("Exact solution with $$a= %1.2f$$", op_c)];
fig_opts = {'Position',  [10, 50, 520, 450]};
figure(fig_opts{:})
fig = gca();
plot_opts = {'Marker','none', 'MarkerSize',30, 'LineStyle','-', 'LineWidth',4};
for i = find(alpha <= 1)
  l_desc = sprintf('$\\alpha = %0.2f$',alpha(i));
  if any(isnan(Lu_ex(i,:,:)))
    u_ex = @(t,x) sol_appr_fcp_cauchy_exact(alpha(i),t,x);
    [Lu_ex(i,:,:),fig] = plot_sol_vs_t(u_ex,x_plot,l_desc,fig,plot_opts,t,x);
  else
    u_exa = @(t,x) squeeze(Lu_ex(i,1:numel(x),1:numel(t)));
    [~,fig] = plot_sol_vs_t(u_exa,x_plot,l_desc,fig,plot_opts,t,x);
  end
end
for i = find(alpha > 1)
  l_desc = sprintf('$\\alpha = %0.2f$',alpha(i));
  u_ex = @(t,x) sol_appr_fcp_cauchy_exact(alpha(i),t,x);
  if any(isnan(Lu_ex(i,:,:)))
    [Lu_ex(i,:,:),fig] = plot_sol_vs_t(u_ex,x_plot,l_desc,fig,plot_opts,t,x);
  else
    u_exa = @(t,x) squeeze(Lu_ex(i,1:numel(x),1:numel(t)));
    [~,fig] = plot_sol_vs_t(u_exa,x_plot,l_desc,fig,plot_opts,t,x);
  end
end
h = legend(fig,'Location','southeast');
set(h,'fontsize',42)
set(gca,'FontSize',39)
xlabel('');
title(desc,'interpreter','latex')

%% Plot approximations
l_desc = @(alpha)sprintf('$\\alpha = %0.2f$',alpha);

for j=1:numel(N)
  desc =[sprintf("Approximate solution with $$N=%d$$:", N(j))];
  fig_opts = {'Position',  [10, 50, 520, 450]};
  figure(fig_opts{:})
  fig = gca();
  for i=find(alpha<=1)
    if any(isnan(Lu_appr(j,i,:,:)))
      df = @(t) df_f(t,x,alpha(i));
      u_appr = @(t,x) sol_appr_fcp_cauchy_R_FD_Dirichlet(alpha(i),N(j),t,x);
      [Lu_appr(j,i,:,:),fig] = plot_sol_vs_t(u_appr,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    else
      u_appra = @(t,x) squeeze(Lu_appr(j,i,1:numel(x),1:numel(t)));
      [~,fig] = plot_sol_vs_t(u_appra,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    end
  end
  for i=find(alpha>1)
    if any(isnan(Lu_appr(j,i,:,:)))
      df = @(t) df_f(t,x,alpha(i));
      u_appr = @(t,x) sol_appr_fcp_cauchy_R_FD_Dirichlet(alpha(i),N(j),t,x);
      [Lu_appr(j,i,:,:),fig] = plot_sol_vs_t(u_appr,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    else
      u_appra = @(t,x) squeeze(Lu_appr(j,i,1:numel(x),1:numel(t)));
      [~,fig] = plot_sol_vs_t(u_appra,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    end
  end
  h = legend(fig);
  set(h,'fontsize',42)
  set(gca,'FontSize',39)
  xlabel('')
  title(desc,'interpreter','latex')
end
save(filename);

%% Plot errors of approximation
l_desc = @(alpha)sprintf('$\\alpha = %0.2f$',alpha);
desc =[];
fprintf('Accuracy of solution at t = 0, x = %1.2f:\n', x_plot);
fprintf('--------------------------------------------------------------------\n');
fprintf('N\\alpha |');
fprintf('%1.2f       |',alpha(alpha<=1),alpha(alpha>1));
fprintf('\n')
fprintf('--------------------------------------------------------------------|\n');
plot_xidx = helper_get_common_grid_idx(x,x,x_plot);
for j=1:numel(N)
  desc =[sprintf("Error of approximate solution with $$N=%d$$:", N(j))];
  fig_opts = {'Position',  [10, 50, 520, 450]};
  figure(fig_opts{:})
  fig = gca();
  fprintf('%8d|',N(j));
  for i=find(alpha<=1)
    u_err = @(t,x) abs(squeeze(Lu_ex(i,1:numel(x),1:numel(t))) - squeeze(Lu_appr(j,i,1:numel(x),1:numel(t))));
    [err,fig] = plot_sol_vs_t(u_err,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    title(desc)
    fprintf('%0.5e|', err(plot_xidx,1));
  end
  for i=find(alpha>1)
    u_err = @(t,x) abs(squeeze(Lu_ex(i,1:numel(x),1:numel(t))) - squeeze(Lu_appr(j,i,1:numel(x),1:numel(t))));
    [err,fig] = plot_sol_vs_t(u_err,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    fprintf('%0.5e|', err(plot_xidx,1));
  end
  set(gca, 'YScale', 'log')
  h = legend(fig);
  set(h,'fontsize',42)
  set(gca,'FontSize',39)
  xlabel('')
  ylim([1e-12,inf])
  title(desc,'interpreter','latex')
  fprintf('\n');
end
fprintf('--------------------------------------------------------------------\n');

%% Local functions
  function [u_sol] = sol_appr_fcp_cauchy_exact(fcp_alpha,t,x)
    u_sol = u_exact_f(t,x);
  end

  function [u_appr] = sol_appr_fcp_cauchy_R_FD_Dirichlet(fcp_alpha,fcp_N,t,x)
    a = max(x);
    fcp_N_hom = fcp_N;        % Discretization parameter for the operator exponential evaluation
    fcp_N_inhom = fcp_N;      % Discretization parameter for the operator exponential evaluation

    %% Fractional Cauchy problem (FCP) propagator parameters
    % Two methods are available:
    % 'hyperbola' - classical hyperbolic contour (viable for t < 10)
    % 'half-hyperbola' - subdivided hyperbolic contour (viable for t < 500)
    % 'auto' - autoselection of the method based on maximum t (this is default)
    method = 'hyperbola';
    fcp_spectral_estimator = 'manual';

    %% Operator Spectrum shape parameters
    fcp_contour_vertex   = NaN;   % NaN - Auto
    fcp_spectral_vertex  = 0;     % Overridden by zero for 'hyperbola' method

    fcp_spectral_angle   = pi/60;
    fcp_critical_vertex  = -pi/6;
    fcp_critical_angle   = pi/2;

    fcp_correction_order = 1;     % Only integer values 0,1 are implemented
    fcp_correction_point = 0;     % NaN - Auto

    fcp_hom = struct('frac_order',fcp_alpha, 'N',fcp_N_hom, ...
      'final_time',T, ...
      'spectral_estimator',fcp_spectral_estimator, ...
      'contour_vertex',fcp_contour_vertex, ...
      'spectral_vertex',fcp_spectral_vertex, ...
      'spectral_angle',fcp_spectral_angle, ...
      'critical_vertex',fcp_critical_vertex, ...
      'critical_angle',fcp_critical_angle, ...
      'correction_order',fcp_correction_order, ...
      'correction_point',fcp_correction_point);
    fcp_inhom = fcp_hom;
    fcp_inhom.N = fcp_N_inhom;

    %% Evaluation of solutions and comparisons
    [u_appr] = fcp_cauchy_hom_appr(t,fcp_hom,u0_res,u1_res,u0,u1);
    [u_appr] = u_appr + fcp_cauchy_inhom_appr(t,fcp_inhom,f0_res,df_res,f0,df,df_strip_height,method);

  end
function r = df_f(t,x,alpha)
  r = 12 * t .* x - 4 * t + 0.2e1 * (-2 + alpha)* (t .^ (1 - alpha))  .* (x .^ 2) .* (x - 1) / gamma((3 - alpha));
end

end
