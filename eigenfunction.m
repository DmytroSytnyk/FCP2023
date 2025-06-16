function [u,du] = eigenfunction(a,k,x)
%
% Copyright (C) 2019-2025, Dmytro Sytnyk. All rights reserved.

%% Function body
  u = sin(pi/a*k .* x);
  if (nargout > 1)
    du = pi/a*k .* cos(pi/a*k .* x);
  end
end


