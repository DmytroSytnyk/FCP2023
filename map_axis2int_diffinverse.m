function [x] = map_axis2int_diffinverse(a,b,y)

%MAP_AXIS2INT_DIFFINVERSE calculates the derivative of map_axis2int_inverse
%
%INPUT:
% a,b - endpoints of the interval
% y   - vector of points from [a, b]
%
%OUTPUT:
% x   - derivatives of log((x - a)./(b - x))
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
x = (b-a)./((y-a).*(b-y));
end

