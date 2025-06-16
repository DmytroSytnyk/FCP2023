function [y] = map_axis2int_diffdirect(a,b,x)

%MAP_AXIS2INT_DIFFDIRECT Calculates the derivative of map_axis2int_direct
%
%INPUT:
% a,b - endpoints of the interval
% x - vector of points from real line
%
%OUTPUT:
% y - vector of derivatives of a + (b-a)*exp(x)./(1 + exp(x))
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
y = (b - a).*exp(x)./((1.0 + exp(x)).^2);
end

