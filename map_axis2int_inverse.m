function [y] = map_axis2int_inverse(a,b,x)

%MAP_AXIS2INT_INVERSE Maps points from the finite interval to real line
%
%INPUT:
% a,b - endpoints of the interval
% x   - vector of points from [a, b]
%
%
%OUTPUT:
% y   - images of the points x from [a,b]
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Function body
% DIGITS = 1024;
Ix = ones(size(x));
% y = vpa(zeros(size(x)),DIGITS);
y = zeros(size(x));
% y = log((x- a*Ix)./(b*Ix-x));
f = (x < (a+b)/2);
% y(f) = vpa(log((x(f)- a*Ix(f))./(b*Ix(f)-x(f))),DIGITS);
y(f) = log((x(f)- a*Ix(f))./(b*Ix(f)-x(f)));
f = (x > (a+b)/2);
% y(f) = vpa(-log((b*Ix(f)-x(f))./(x(f)-a*Ix(f))),DIGITS);
y(f) = -log((b*Ix(f)-x(f))./(x(f)-a*Ix(f)));
% f = (isnan(y) | isinf(y));
% if any(f)
%   y(f) = -log((b*ones(nnz(f),1) - x(f))./(x(f) - a*ones(nnz(f),1)));
% end
y = double(y);
end

