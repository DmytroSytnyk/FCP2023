function [err,fig] = plot_err_vs_n(fsol_err,alpha,N,s_desc,fig,plot_opts,t,x,err_norm_t,err_norm_x)

%PLOT_ERR_VS_N Plot the error norm of solution with respect to N
%
%INPUT:
% fsol_err   - Function handle err = @(alpha,N,t,x)
% N          - Discretization parameter for the operator exponential evaluation
% alpha      - Fractional derivative order parameter
% s_desc     - Additional description for the figure  
% fig        - Figure to plot on (empty means generate new figure)
% plot_opts  - Options for the plot
% t          - Grid in time 
% x          - Grid in space
% err_norm_t - Error norm in time
% err_norm_x - Error norm in space
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin==4) 
  fig = [];
  plot_opts = [];
elseif (nargin==5) 
  plot_opts = [];
end
if (nargin <= 6)
  % Grids in space and time
  t = [0:1/100:1];
  x = [];
elseif (nargin == 7)
  x = [];
end  
if (nargin <= 8)
  % Error norm 
  err_norm_t = Inf;
  err_norm_x = Inf;
end
if (nargin == 9)
  err_norm_x = Inf;
end

if isempty(plot_opts)
  plot_opts = {};
end

if isempty(x)
  x = (0:0.01:1).';
end

if isempty(fig)
  figure;
  fig = gca;
else
  axes(fig);
end

%% Function body
parameters = {N,alpha};
combs = combvec(parameters{:}).'; % Use cells as arguments
nr = size(combs,1);

% Main calculation cycle
p = zeros(nr,1);
err = zeros(numel(alpha),numel(N));
for i =1:nr
  fcp_N = combs(i,1);      % Discretization parameter for the operator exponential evaluation
  fcp_alpha = combs(i,2);  % Fractional derivative order parameter
  % Solve problem and obtain error
  u_err = fsol_err(fcp_alpha,fcp_N,t,x);
  err(find(alpha==fcp_alpha),find(N==fcp_N))= vecnorm(vecnorm(u_err,err_norm_x,1),err_norm_t);
end

% Plot the results
s_info = '';
set(gca, 'YScale', 'log')
for i = 1:numel(alpha)
  semilogy(N,err(i,:),'DisplayName',sprintf("$\\alpha$=%0.2f",alpha(i)),plot_opts{:});
  hold on
end
hold off
grid on
title({'Accuracy of solution vs N';s_desc;s_info});
hl = legend('show');
set(hl, 'Interpreter','latex')
xlabel('N'); ylabel('Error');

end
