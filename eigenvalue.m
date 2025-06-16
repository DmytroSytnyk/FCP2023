function [ev] = eigenvalue(a,k, op_c)
%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved.

%% Parse input parameters
if (nargin == 2)
  op_c = 1;
elseif (nargin == 3)
  % Do nothing
else
  error('Wrong number of input arguments.');
end

%% Function body
ev =  op_c*(pi/a) ^ 2 * k .^ 2;
end
