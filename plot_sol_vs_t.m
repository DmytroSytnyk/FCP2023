function [u_sol, fig] = plot_sol_vs_t(fsol,x_plot,l_desc,fig,plot_opts,t,x)

%PLOT_SOL_VS_T Plot solution with respect to time
%
%INPUT:
% fsol_err   - function handle err = @(t,x)
% l_desc     - Description for the legend
% fig        - Figure to plot on (empty means generate new figure)
% plot_opts  - Options for the plot
% t          - Grid in time
% x          - Grid in space
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin<=3)
  fig = [];
elseif (nargin<=4)
  plot_opts = [];
end
if (nargin <= 5)
  % Grids in space and time
  t = [0:1/100:1];
  x = [];
elseif (nargin == 6)
  x = [];
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
u_sol = fsol(t,x);

% Plot the results
plot_xidx = helper_get_common_grid_idx(x,x,x_plot);
s_info = '';
for i = numel(plot_xidx)
%   plot(t,real(u_sol(plot_xidx(i),:)),'DisplayName',sprintf("%s(t,%0.1f)",l_desc,x(plot_xidx(i))),plot_opts{:});
  plot(fig,t,real(u_sol(plot_xidx(i),:)),'DisplayName',l_desc,plot_opts{:});
  hold on
end
grid on
% title({'Exact vs approximate solution over time';s_desc;s_info});
hl = legend('show');
set(hl, 'Interpreter','latex')
xlabel('t'); ylabel('');

end
