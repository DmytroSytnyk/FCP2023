function [r] = RLI_approx(f, alpha, t, a, b, N, d, tune_N2t)

%RLI_APPROX approximates the Riemann-Liouville integral of the function defined on [a,b] using Sinc Quadrature rule
%
% [r] = RLI_approx(f, alpha, t, a, b, N, d)
% approximates the fractional Riemann-Liouville integral
%    r  = 1/Gamma(alpha)*int((t-x)^alpha*f(s),s=a..t)
% using the Sinc-quadrature formula composed with the conformal mapping of
% (a,t) into (-Inf, Inf)
%
%INPUT:
% f        - a vectorizable function y = @(x) that return values of the 
%            integrand. It is assumed that f is vectorizable and for each 
%            argument returns a column vector.
% alpha    - fractional order parameter
% t        - values of integral's upper bound (can be an array)
% a        - left endpoint of the inteval
% b        - right endpoint of the interval
% N        - discretization parameter
% d        - half-height of the analyticity domain for f(psi(x))
% tune_N2t - automatically tune N to match accuracy for different t
%
%OTUPUT:
% r        - results of approximation vector of size = seize(f(t))
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input arguments

if (nargin == 7)
  tune_N2t = false;
elseif (nargin == 8)
  % Do nothing
else
  error('Wrong number of input arguments.');
end

if (d<=0) || (N<0)
  error('Input parameters are out of range (d>0, N > 0)');
end
if ~(a < b) || (any(t > b | t < a))
  error('Ensure b > a and all values t from interval [a, b]');
end

%% Function body
% The next function defines the part of the integrand that comes from the
% product of 1/(t-s)^alpha by the derivative of mapping s = a
% +t*exp(x)/(1+exp(x)) with the term t^(alpha+1) factored out in front of
% the integral
% sk1 = @(x) exp(x)./((1+exp(x)).^(alpha+1));
sk2 = @(x) (exp(-x/(alpha+1))+exp(x*alpha/(alpha+1))).^(-(alpha+1));
% Evaluate the RL integral for t(1) and determine the size of output r
%The decay of the integrand is ~e^(-x), ~e(-alpha*x) hence the alpha in the
%estimate of h (see Eq. (3.2) from F. Stenger "Summary of sinc methods")
g = gamma(alpha);
eps = min(1, alpha);
t0 = max(t);
nt = numel(t);
[X,x,h] = gen_grid(nt);
% S = (sk1(x)).';
S = (sk2(x)).';
% Handle the case when f is a constant
fX1 = f(X(1));
fX12 = f(X(1:2));
if (size(fX1) == size(fX12)) 
  vf = repmat(fX1,numel(X));
  rl = (h*t(nt)^alpha/g)*vf*S;
  r = zeros([size(rl),numel(t)]);
  r(:,nt) = rl;
  for i = 1:(nt-1)
    if (tune_N2t)
      [X,x] = gen_local_grid(i);
      % S = (sk1(x)).';
      S = (sk2(x)).';
    end
    if (t(i) ~= a)
      r(:,i) = (h*t(i)^alpha/g)*vf(:,numel(X))*S;
    end
  end
else
  vf = f(X);
  rl = (h*t(nt)^alpha/g)*vf*S;
  r = zeros([size(rl),nt]);
  r(:,nt) = rl;
  for i=1:(nt-1)
    if (tune_N2t)
      [X,x,h] = gen_grid(i);
      % S = (sk1(x)).';
      S = (sk2(x)).';
    else
      X = map_axis2int_direct(a,t(i),x);
    end
    r(:,i) = (h*t(i)^alpha/g)*f(X)*S;
  end
end
r = reshape(squeeze(r),size(repmat(fX1,size(t))));

function [N] =  RLI_adj_N2t(s,N0,t0,t) 
%   N = round((sqrt(N0) + s.*log(t/t0)).^2);
  N = round((max(1,sqrt(N0) + s.*log(t/t0))).^2);
%   N = N0;
%   N = round((max(4,sqrt(N0) + s.*(log(t)-log(t0)))).^2);
  if (mod(N0,2) ~= mod(N,2))
    N = N +1;
  end
end

function [X,x,h] = gen_grid(k)
  % Adjust N so that the accuracy whould be approximatelly equal for all t
  Nt = RLI_adj_N2t(alpha/sqrt(2*pi*d*eps),N,t0,t(k));
  N1 = round(eps*Nt);
  N2 = round(min(1/alpha, 1)*Nt);
  h = quad_SQ_get_h(Nt,d,eps);
  x = h*(-N1:N2);
  X = map_axis2int_direct(a,t(k),x);
end
end
