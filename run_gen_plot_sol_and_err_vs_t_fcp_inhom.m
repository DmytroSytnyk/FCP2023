function [] = run_gen_plot_sol_and_err_vs_t_fcp_inhom()
%
% Copyright (C) 2022 - 2025, Dmytro Sytnyk. All rights reserved. 

environment

% Lists of parameters to perform the calculations on
% The following code will be executed for each combination of the
% parameters from the list
% alpha = [0.1, 0.5, 0.7, 0.9, 1, 1.1, 1.5, 1.7, 1.9];
%  alpha = [1, 0.7, 0.5, 0.3, 0.1];
alpha = [1, 1.2, 1.5, 1.7, 1.9];

% N = [448];
N = [32,64,128];
% N = [8,16,32,64,128,256,512,1024];

%% Grid parameters
T = 5;
Nt = 200;
Nx = 100;
t = [0:T/Nt:T];
% t = [0];
x = (0:1/Nx:1).';
x_plot = 0.5;

%% Exact solution parameters
ex_NI = 256;

%% Problem parameters
%NOTE!
% For the purpose of this test both f0 and df should be constant
ev_f0 = 1;                       % Eigenfunction index for f0
f0_c = 1;                        % Value of the right-hand side at zero
ev_df = 4;                       % Eigenfunction index for df
df_c = 1;
df_t = @(t) df_c*ones(size(t));  % Derivative of the right-hand side
% This is set to maximum pi/4 since df_t is constant
df_strip_height = pi/4;          % The eye-shaped domain analyticity angle

%% Resolvent parameters
op_c = 1;                        % We consider A = -op_c*d^2/dx^2

%% Output file name
file_prefix = 'run_gen_plot_sol_and_err_vs_t_fcp_inhom_';
salpha = ['alpha_',sprintf('%1.1f_', sort(alpha))];
sN = ['N_',sprintf('%d_', sort(N))];
st = sprintf('t_%1.0f_%1.0f_',min(t),max(t));
sev = sprintf('nf0_%d_ndf_%d_f0c_%1.0f_dfc_%1.0f_',ev_f0,ev_df,f0_c,df_c);
filename = [file_prefix,salpha,sN,st,sev,'.mat'];

%% Load data file
Lu_ex = NaN(numel(alpha),numel(x),numel(t));
Lu_appr = NaN(numel(N),numel(alpha),numel(x),numel(t));
if exist(filename, 'file')
  fprintf("Associated data file:\n %s \n found.\n Loading ... ",filename);
  load(filename)
  fprintf("success.\n");
end

%% Plot exact solutions
desc =[sprintf("Exact solution with $$a= %1.2f$$, $$N_I=%d$$", op_c, ex_NI)];
fig_opts = {'Position',  [10, 50, 520, 450]};
figure(fig_opts{:})
fig = gca();
plot_opts = {'Marker','none', 'MarkerSize',30, 'LineStyle','-', 'LineWidth',4};
for i = find(alpha <= 1)
  l_desc = sprintf('$\\alpha = %0.2f$',alpha(i));
  if any(isnan(Lu_ex(i,:,:)))
    u_ex = @(t,x) sol_ex_fcp_inhom_R_exact_eigenfunc(alpha(i),t,x);
    [Lu_ex(i,:,:),fig] = plot_sol_vs_t(u_ex,x_plot,l_desc,fig,plot_opts,t,x);
  else
    u_exa = @(t,x) squeeze(Lu_ex(i,1:numel(x),1:numel(t)));
    plot_xidx = helper_get_common_grid_idx(x,x,x_plot);
    [~,fig] = plot_sol_vs_t(u_exa,x_plot,l_desc,fig,plot_opts,t,x);
  end
end
for i = find(alpha >= 1)
  l_desc = sprintf('$\\alpha = %0.2f$',alpha(i));
  u_ex = @(t,x) sol_ex_fcp_inhom_R_exact_eigenfunc(alpha(i),t,x);
  if any(isnan(Lu_ex(i,:,:)))
    [Lu_ex(i,:,:),fig] = plot_sol_vs_t(u_ex,x_plot,l_desc,fig,plot_opts,t,x);
  else
    u_exa = @(t,x) squeeze(Lu_ex(i,1:numel(x),1:numel(t)));
    plot_xidx = helper_get_common_grid_idx(x,x,x_plot);
    [~,fig] = plot_sol_vs_t(u_exa,x_plot,l_desc,fig,plot_opts,t,x);
  end
end
h = legend(fig,'Location','southeast');
set(h,'fontsize',42)
set(gca,'FontSize',39)
% title("Exact solution")
title(desc,'interpreter','latex')

%% Plot approximations
l_desc = @(alpha)sprintf('$\\alpha = %0.2f$',alpha);

for j=1:numel(N)
  desc =[sprintf("Approximate solution with $$N=%d$$:", N(j)); sprintf("$$f_0 = %1.2f \\sin{(%d \\pi x)}$$, $$f'(t) = %1.2f \\sin{(%d \\pi x)}$$",f0_c,ev_f0,df_c,ev_df)];
  fig_opts = {'Position',  [10, 50, 520, 450]};
  figure(fig_opts{:})
  fig = gca();
  for i=find(alpha<=1)
    if any(isnan(Lu_appr(j,i,:,:)))
      u_appr = @(t,x) sol_appr_fcp_inhom_R_exact_eigenfunc(alpha(i),N(j),t,x);
      [Lu_appr(j,i,:,:),fig] = plot_sol_vs_t(u_appr,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    else
      u_appra = @(t,x) squeeze(Lu_appr(j,i,1:numel(x),1:numel(t)));
      [~,fig] = plot_sol_vs_t(u_appra,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    end

  end
  for i=find(alpha>1)
    if any(isnan(Lu_appr(j,i,:,:)))
      u_appr = @(t,x) sol_appr_fcp_inhom_R_exact_eigenfunc(alpha(i),N(j),t,x);
      [Lu_appr(j,i,:,:),fig] = plot_sol_vs_t(u_appr,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    else
      u_appra = @(t,x) squeeze(Lu_appr(j,i,1:numel(x),1:numel(t)));
      [~,fig] = plot_sol_vs_t(u_appra,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    end
  end
  h = legend(fig);
  set(h,'fontsize',42)
  set(gca,'FontSize',39)
  title(desc,'interpreter','latex')
end
save(filename);

%% Plot errors of approximation
l_desc = @(alpha)sprintf('$\\alpha = %0.2f$',alpha);
for j=1:numel(N)
  desc =[sprintf("Error of the approximate solution with $$N=%d$$:", N(j)); sprintf("$$f_0 = %1.2f \\sin{(%d \\pi x)}$$, $$f'(t) = %1.2f \\sin{(%d \\pi x)}$$",f0_c,ev_f0,df_c,ev_df)];
  fig_opts = {'Position',  [10, 50, 520, 450]};
  figure(fig_opts{:})
  fig = gca();
  for i=find(alpha<=1)
    u_err = @(t,x) abs(squeeze(Lu_ex(i,1:numel(x),1:numel(t))) - squeeze(Lu_appr(j,i,1:numel(x),1:numel(t))));
    [~,fig] = plot_sol_vs_t(u_err,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
    title(desc,'interpreter','latex')
  end
  for i=find(alpha>1)
    u_err = @(t,x) abs(squeeze(Lu_ex(i,1:numel(x),1:numel(t))) - squeeze(Lu_appr(j,i,1:numel(x),1:numel(t))));
    [~,fig] = plot_sol_vs_t(u_err,x_plot,l_desc(alpha(i)),fig,plot_opts,t,x);
  end
  set(gca, 'YScale', 'log')
  set(h,'fontsize',42)
  set(gca,'FontSize',39)

end

  function [u_sol] = sol_ex_fcp_inhom_R_exact_eigenfunc(fcp_alpha,t,x)
    a = max(x);
    ex_NJ = ex_NI/min(1,fcp_alpha);
    u_sol = fcp_cauchy_inhom_ex_sol_eigenfunction(t, fcp_alpha, x, op_c, ev_f0, f0_c, ev_df, df_t, ex_NI, ex_NJ, df_strip_height);
  end

  function [u_appr] = sol_appr_fcp_inhom_R_exact_eigenfunc(fcp_alpha,fcp_N,t,x)
    method = 'hyperbola';
    fcp_spectral_estimator = 'manual';

    % Operator Spectrum shape parameters
    fcp_contour_vertex   = NaN;   % NaN - Auto
    fcp_spectral_vertex  = 0;     % Overridden by zero for 'hyperbola' method
    fcp_spectral_angle  = pi/60;
    fcp_critical_vertex = -pi/6;
    fcp_critical_angle  = pi/2;
    fcp_correction_order = 1;     % Only integer values 0,1 are implemented
    fcp_correction_point = 0;     % NaN - Auto
    fcp = struct('frac_order',fcp_alpha, 'N',fcp_N, ...
      'final_time',T, ...
      'spectral_estimator',fcp_spectral_estimator, ...
      'contour_vertex',fcp_contour_vertex, ...
      'spectral_vertex',fcp_spectral_vertex, ...
      'spectral_angle',fcp_spectral_angle, ...
      'critical_vertex',fcp_critical_vertex, ...
      'critical_angle',fcp_critical_angle, ...
      'correction_order',fcp_correction_order, ...
      'correction_point',fcp_correction_point);

    % Resolvent of spatial operator and the right hand side of the problem
    a = max(x);
    f0_res = @(z,u) resolvent_exact_eigenfunc(z,ev_f0,x,a,op_c,u);
    df_res = @(z,u) resolvent_exact_eigenfunc(z,ev_df,x,a,op_c,u);
    f0 = f0_c*eigenfunction (a, ev_f0, x);
    df = @(t) df_t(t).*eigenfunction (a, ev_df, x);

    [u_appr] = fcp_cauchy_inhom_appr(t,fcp,f0_res,df_res,f0,df,df_strip_height,method);
  end
end