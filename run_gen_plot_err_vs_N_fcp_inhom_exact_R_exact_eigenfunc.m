function [alpha, N, err] = run_gen_plot_err_vs_N_fcp_inhom_exact_R_exact_eigenfunc()
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

environment

%% Error norms
err_norm_x = Inf;
err_norm_t = Inf;

% Lists of parameters to perform the calculations on
% The following code will be executed for each combination of the
% parameters from the list
% alpha = [0.1, 0.5, 0.7, 0.9, 1, 1.1, 1.5, 1.7, 1.9];
alpha = [1, 0.7, 1.2, 0.5, 1.5, 0.3, 1.7, 0.1, 1.9];
% alpha = [0.1, 0.4, 1, 1.4, 1.9];
%  alpha = [1, 0.7, 0.5, 0.3, 0.1];
N = [32:32:2048];
% N = power(2,5:9);

%% Grid parameters
Nt = 200;
Nx = 100;
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
op_c = 1;               % We consider A = -op_c*d^2u/dx^2

%% Output file name
file_prefix = 'run_gen_plot_err_vs_N_fcp_inhom_exact_R_exact_eigenfunc_';

%% Script body
method = 'hyperbola';
alpha = unique(alpha,'stable');
% desc=sprintf("D= %.2f, u_0,u_1 - %d,%d-th eigenfunction.\n Contour: %s",op_c, ev_n0,ev_n1,method);
filename="";
storage=struct('N',[]);
desc = [];
fig_opts = {'Position',  [10, 50, 520, 450]};
figure(fig_opts{:})
fig = gca();

T = 1; t = [0,T];
%% Load data file
Lu_ex = NaN(numel(alpha),numel(x),numel(t));
Lu_appr = NaN(numel(N),numel(alpha),numel(x),numel(t));
load_data();

plot_opts = {'Marker','.', 'MarkerSize',30, 'LineStyle','-', 'LineWidth',3};
[~,fig] = plot_err_vs_n(@error_fcp_inhom_exact_R_exact_eigenfunc,alpha(alpha<=1),N,desc,fig,plot_opts,t,x,err_norm_t,err_norm_x);
xlim([0,max(N)+min(N)])
ylim([1e-16,1])
h = legend(fig);
set(h,'fontsize',32)
set(gca,'FontSize',26)
xlabel(''); ylabel('');
% hold on
save_data();
saveas(gcf, [filename,'.eps']);

T = 5; t = [0,T];
% figure
hold on
%% Load data file
Lu_ex = NaN(numel(alpha),numel(x),numel(t));
Lu_appr = NaN(numel(N),numel(alpha),numel(x),numel(t));
load_data();
plot_opts = {'Marker','x', 'MarkerSize',11, 'LineStyle','-', 'LineWidth',3};
[~,fig] = plot_err_vs_n(@error_fcp_inhom_exact_R_exact_eigenfunc,alpha(alpha>=1),N,desc,fig,plot_opts,t,x,err_norm_t,err_norm_x);

xlim([0,max(N)+min(N)])
ylim([1e-16,1])
h = legend(fig);
set(h,'fontsize',32)
set(gca,'FontSize',26)
xlabel(''); ylabel('');
save_data();
saveas(gcf, [filename,'.eps']);

%% Local functions
  % Calculates error of exact and approximate solution using stored data,
  % whenever avaliable 
  function [err] = error_fcp_inhom_exact_R_exact_eigenfunc(fcp_alpha,fcp_N,t,x)
    alpha_idx = find(alpha == fcp_alpha);
    N_idx = find(N == fcp_N);
    if any(isnan(Lu_ex(alpha_idx,:,:)))
      Lu_ex(alpha_idx,:,:) = sol_ex_fcp_inhom_R_exact_eigenfunc(alpha(alpha_idx),t,x);
    else
      fprintf('Exact solution data for alpha = %1.1f found. Skipping.\n',fcp_alpha);
    end
    if any(isnan(Lu_appr(N_idx,alpha_idx,:,:)))
      fprintf('Approximate solution data for alpha = %1.1f, N = %d NOT found. Calculating.\n',fcp_alpha,fcp_N);
      Lu_appr(N_idx,alpha_idx,:,:) =  sol_appr_fcp_inhom_R_exact_eigenfunc(alpha(alpha_idx),N(N_idx),t,x);
      save_data();
    else
      fprintf('Approximate solution data for alpha = %1.1f, N = %d found. Skipping.\n',fcp_alpha,fcp_N);
    end
    err = squeeze(Lu_ex(alpha_idx,:,:)) - squeeze(Lu_appr(N_idx,alpha_idx,:,:));
  end

  % Evaluates exact solution
  function [u_sol] = sol_ex_fcp_inhom_R_exact_eigenfunc(fcp_alpha,t,x)
    a = max(x);
    ex_NJ = ex_NI/min(1,fcp_alpha);
    u_sol = fcp_cauchy_inhom_ex_sol_eigenfunction(t, fcp_alpha, x, op_c, ev_f0, f0_c, ev_df, df_t, ex_NI, ex_NJ, df_strip_height);
  end

  % Calculates approximate solution
  function [u_appr] = sol_appr_fcp_inhom_R_exact_eigenfunc(fcp_alpha,fcp_N,t,x)
    method = 'hyperbola';
    fcp_spectral_estimator = 'manual';

    % Operator Spectrum shape parameters
    fcp_contour_vertex   = NaN;   % NaN - Auto
    fcp_spectral_vertex  = 0;     % Overridden by zero for 'hyperbola' method
    fcp_spectral_angle   = pi/60;
    fcp_critical_vertex  = -pi/6;
    fcp_critical_angle   = pi/2;
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

  % Initialize storage using file (if present)
  function  load_data()
    salpha = ['alpha_',sprintf('%1.1f_', sort(alpha))];
    sN = ['N_',sprintf('%d_', [min(N),max(N)])];
    st = sprintf('t_%1.0f_%1.0f_',min(t),max(t));
    sev = sprintf('nf0_%d_ndf_%d',ev_f0,ev_df);
    filename = [file_prefix,salpha,sN,st,sev];
    fprintf('Attemption to load data from file:\n  %s\n',filename);
    if isfile([filename,'.mat'])
      fprintf("Associated data file found. Loading ... ");
      storage = load([filename,'.mat'],'Lu_ex','Lu_appr','N','alpha','x','t','ex_NI');
      copyfile([filename,'.mat'],[filename,'.bak'])
      [mN,mNi] = ismember(N,storage.N);
      if any(mN) && isequal(x,storage.x) && isequal(t,storage.t)
        Lu_ex = storage.Lu_ex;
        Lu_appr(mN,:,:,:) = storage.Lu_appr(mNi(mN),:,:,:);
        fprintf("succes.\n");
      else
        fprintf("Provided N,x,t are incompatible with stored data.")
      end
    else
      storage = struct('N',N,'t',[min(t),max(t)],'x',x,'alpha',alpha,'ex_NI',ex_NI);
    end
  end

  % Save data to storage and file
  function save_data()
    [mN,mNi] = ismember(N, storage.N);
    if any(mN) && isequal(x,storage.x) && isequal(t,storage.t)
      storage.Lu_ex = Lu_ex;
      storage.Lu_appr(mNi(mN),:,:,:) = Lu_appr(mN,:,:,:);
      save([filename,'.mat'],'-struct','storage');
      fprintf("succes.\n");
    else
      fprintf("Provided N,x,t are incompatible with stored data.")
    end
  end
end
