function [alpha, N, err] = run_gen_plot_err_vs_N_fcp_cauchy_R_FD_Dirichlet()
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

environment

%% Error norms
err_norm_x = Inf;
err_norm_t = Inf;

% Lists of parameters to perform the calculations on
% The following conde will be executed for each combination of the
% parameters from the list
% alpha = [0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.5, 1.7, 1.9];
alpha = [1, 0.7, 1.2, 0.5, 1.5, 0.3, 1.7, 0.1, 1.9];
% alpha = [0.1, 0.4, 1, 1.4, 1.9];
%  alpha = [1, 0.7, 0.5, 0.3, 0.1];
N = [32:32:512];
% N = power(2,5:6);


%% Grid parameters
T = 1;
Nt = 200;
Nx = 1e4;
t = [0:T/Nt:T];
% t = [0];
x = (0:1/Nx:1).';

%% Exact solution parameters
% u_exact_f = @(t,x) x .* (x - t .^ 2 + 0.1e1 / 0.2e1) .* (x - 0.1e1);
u_exact_f = @(t,x) x.^2 .* (x - t .^ 2 + 0.1e1 / 0.2e1) .* (x - 0.1e1);

%% Problem parameters
% Operator parameters
op_c = 1;             % We consider A = -op_c*d^2u/dx^2

% Initial values
u0 = u_exact_f(0,x);
u1 = zeros(size(x));

% Right-hand side
% Below df(x) means f(0,x) and df_f(t,x,alpha) means f'(t,x)
f0_f = @(x) -12 * x .^ 2 + 3 * x + 1;
f0 = f0_f(x);
df_strip_height = pi/4; % The eye-shaped domain analyticity angle for df
  
%% Resolvent of spatial operator and the right hand side of the problem
a = max(x);
u0_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);
u1_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);
f0_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);
df_res = @(z,u) resolvent_fd_1d_epdo_bc_dirichlet(z,x,u,op_c,0,0,[0,a]);


fprintf('|A u_0| = %1.3e , |A df| = %1.3e .\n',vecnorm(u0_res(0,u0),err_norm_x,1),vecnorm(df_res(0,df_f(1,x,0.5)),err_norm_x,1));

%% Output file name
file_prefix = 'run_gen_plot_err_vs_N_fcp_cauchy_R_FD_Dirichlet_';

%% Script body
method = 'hyperbola';
alpha = unique(alpha,'stable');
% desc=sprintf("D= %.2f, u_0,u_1 - %d,%d-th eigenfunction.\n Contour: %s",op_c, ev_n0,ev_n1,method);
filename="";
storage=struct('N',[]);
desc = [];
fig_opts = {'Position',[10, 50, 520, 450]};
figure(fig_opts{:})
fig = gca();

T = 1; t = [0,T];
%% Load data file
Lu_ex = NaN(numel(alpha),numel(x),numel(t));
Lu_appr = NaN(numel(N),numel(alpha),numel(x),numel(t));
load_data();

plot_opts = {'Marker','.', 'MarkerSize',30, 'LineStyle','-', 'LineWidth',3};
[~,fig] = plot_err_vs_n(@error_appr_fcp_cauchy_R_FD_Dirichlet,alpha(alpha<=1),N,desc,fig,plot_opts,t,x,err_norm_t,err_norm_x);
xlim([0,max(N)+min(N)])
% ylim([1e-16,1])
h = legend(fig);
set(h,'fontsize',32)
set(gca,'FontSize',26)
xlabel(''); ylabel('');
% hold on
save_data();
saveas(gcf,[filename,'.eps']);

T = 1;t = [0,T];
% figure
hold on
%% Load data file
Lu_ex = NaN(numel(alpha),numel(x),numel(t));
Lu_appr = NaN(numel(N),numel(alpha),numel(x),numel(t));
load_data();
plot_opts = {'Marker','x', 'MarkerSize',11, 'LineStyle','-', 'LineWidth',3};
[~,fig] = plot_err_vs_n(@error_appr_fcp_cauchy_R_FD_Dirichlet,alpha(alpha>1),N,desc,fig,plot_opts,t,x,err_norm_t,err_norm_x);

xlim([0,max(N)+min(N)])
% ylim([1e-16,1])
h = legend(fig);
set(h,'FontSize',32)
set(gca,'FontSize',26)
xlabel(''); ylabel('');
save_data();
saveas(gcf, [filename,'.eps']);

%% Local functions
  % Calculates error of exact and approximate solution using stored data,
  % whenever avaliable
  function [err] = error_appr_fcp_cauchy_R_FD_Dirichlet(fcp_alpha,fcp_N,t,x)
    alpha_idx = find(alpha == fcp_alpha);
    N_idx = find(N == fcp_N);
    if any(isnan(Lu_ex(alpha_idx,:,:)))
      Lu_ex(alpha_idx,:,:) = sol_ex_fcp_cauchy_R_FD_Dirichlet(alpha(alpha_idx),t,x);
    else
      fprintf('Exact solution data for alpha = %1.1f found. Skipping.\n',fcp_alpha);
    end
    if any(isnan(Lu_appr(N_idx,alpha_idx,:,:)))
      fprintf('Approximate solution data for alpha = %1.1f, N = %d NOT found. Calculating.\n',fcp_alpha,fcp_N);
      Lu_appr(N_idx,alpha_idx,:,:) =  sol_appr_fcp_cauchy_R_FD_Dirichlet(alpha(alpha_idx),N(N_idx),t,x);
      save_data();
    else
      fprintf('Approximate solution data for alpha = %1.1f, N = %d found. Skipping.\n',fcp_alpha,fcp_N);
    end
    err    = squeeze(Lu_ex(alpha_idx,:,:)) - squeeze(Lu_appr(N_idx,alpha_idx,:,:));
  end

  % Evaluates exact solution
  function [u_sol] = sol_ex_fcp_cauchy_R_FD_Dirichlet(fcp_alpha,t,x)
    u_sol = u_exact_f(t,x);
  end

  % Evaluates approximate solution
  function [u_appr] = sol_appr_fcp_cauchy_R_FD_Dirichlet(fcp_alpha,fcp_N,t,x)
    
    df = @(t) df_f(t,x,fcp_alpha); 
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

    fcp_spectral_angle  = pi/60;
    fcp_critical_vertex = -pi/6;
    fcp_critical_angle  = pi/2;

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
    % fcp_cauchy_hom_appr(t,fcp,u0_res,u1_res,u0,u1,method)
    [u_appr] = fcp_cauchy_hom_appr(t,fcp_hom,u0_res,u1_res,u0,u1);
    % fcp_cauchy_inhom_appr(t,fcp,f0_res,df_res,f0,df,df_strip_height, method, smoothness_f0,  smoothness_df)
    [u_appr] = u_appr + fcp_cauchy_inhom_appr(t,fcp_inhom,f0_res,df_res,f0,df,df_strip_height,method);

  end
  function r = df_f(t,x,alpha)
    %   r = (12 * t .* x) - (4 * t) + 0.2e1 * (-2 + alpha) / gamma((3 - alpha)) * (t .^ (1 - alpha)) .* (x .^ 2) .* (x - 1) ;
    r = 12 * t .* x - 4 * t + 0.2e1 * (-2 + alpha)* (t .^ (1 - alpha))  .* (x .^ 2) .* (x - 1) / gamma((3 - alpha));
    %   r = 2* (-2 + alpha) * t .^ (1 - alpha)  .* x .* (x - 1) / gamma((3 - alpha)) + 4 * t;
  end

% Initialize storage using file (if present)
  function  load_data()
    salpha = ['alpha_',sprintf('%1.1f_', sort(alpha))];
    sN = ['N_',sprintf('%d_', sort(N))];
    st = sprintf('t_%1.0f_%1.0f_',min(t),max(t));
    sev = sprintf('Nt_%d_Nx_%d',Nt,Nx);
    filename = [file_prefix,salpha,sN,st,sev];
    fprintf('Attemption to load data from file:\n  %s\n',filename);
    if isfile([filename,'.mat'])
      fprintf("Associated data file found. Loading ... ");
      storage = load([filename,'.mat'],'Lu_ex','Lu_appr','N','alpha','x','t');
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
      storage = struct('N',N,'t',[min(t),max(t)],'x',x,'alpha',alpha);
    end
  end

% Save data to storage and file
  function  save_data()
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
