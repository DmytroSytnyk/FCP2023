function [y] = map_axis2int_direct(a,b,x)

%MAP_AXIS2INT_DIRECT Maps points from the real line to finite interval
%
%INPUT:
% a,b - endpoints of the interval
% x - vector of points from real line
%
%OUTPUT:
% y - vector of images of x that lie on [a, b]
%
% Copyright (C) 2022, Dmytro Sytnyk. All rights reserved. 

%% Function body
% y = a + (b-a)*exp(x)./(1 + exp(x));
y = a + (b-a)./(1 + exp(-x));
% Split version that could give better stability
% m = x>=0;
% y(m)= a + (b-a)./(1 + exp(-x(m)));
% m = x<0;
% y(m) = a + (b-a)*exp(x(m))./(1 + exp(x(m)));
end

