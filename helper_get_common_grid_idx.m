function [g1idx,g2idx,grid_eps_out] = helper_get_common_grid_idx(grid1,grid2,x,grid_eps)

%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin < 4) || isempty(grid_eps)
  grid_eps = NaN;
end
grid1 = reshape(grid1,1,[]);
grid2 = reshape(grid2,1,[]);
x = reshape(x,1,[]);

%% Function body
n = min(1e4,length(x));
i = 1;
g1idx = zeros(1,numel(x));
while i <= max(1,length(x))
  si = min(i+n-1,numel(x));
  A = repmat(grid1,si-i+1,1);
  [~,g1idx(i:si)] = min(abs(A.'- x(i:si)));
%   g1idx(i:si) = g1idx(i:si) + ones(1,si+1-i)*(i-1);
  i = i+n;
end
i = 1;
g2idx = zeros(1,numel(x));
while i <= max(1,length(x))
  si = min(i+n-1,numel(x));
  A = repmat(grid2,si-i+1,1);
  [~,g2idx(i:si)] = min(abs(A.'- x(i:si)));
%   g2idx(i:si) = g2idx(i:si) + ones(1,si+1-i)*(i-1);
  i = i+n;
end
grid_eps_out = max(abs(grid1(g1idx)-grid2(g2idx)));
if ~isnan(grid_eps) && (grid_eps_out>grid_eps)
  error('Provided grids are incompatible: max diff is %.3e',max(abs(grid1(g1idx)-grid2(g2idx))));
end
end
