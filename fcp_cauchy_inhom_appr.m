function [fcp_prop_t, Sf0, Sf0_matr_func] = fcp_cauchy_inhom_appr(t, fcp, f0_res, df_res, f0, df, df_strip_height, method, smoothness_f0,  smoothness_df)

% FCP_CAUCHY_INHOM_APPR approximate solution of inhomogeneous fractional Cauchy problem
%
%INPUT:
% t                 - the array time-points for the solution approximation
% fcp               - fractional Cauchy problem (FCP) structure that 
%                     defines problem  and operator parameters. 
%   Structure fcp is assumed to have the following fields:
%      frac_order      : fractional derivative order parameter
%      spectral_angle  : actual spectral angle
%      spectral_vertex : actual spectral vertex
%      critical_angle  : critical angle
%      critical_vertex : critical vertex
%      contour_vertex  : the prescribed value of contour at zero
%      final_time      : final time for the approximation
% f0_res           - function implementing the evaluation of resolvent 
%                    applied to f0
% df_res           - function implementing the evaluation of resolvent 
%                    applied to df
%   For FCP the resolvent has to be defined as
%                   R(z) = z^(alpha-1)*(z^alpha*I + A)^(-1) 
%   whence the minus sign in the definition of the contour below.
%   So in order to re-use the resolvent functions that implement convenient
%   resolvent implementation  R(z) = (z*I - A)^(-1) we will define
%   adhoc function wrapper. This will allow to freely reuse the resolvent
%   implementations across different problems.
%
% f0               - the value of RHS at 0 f(0)
% df               - the function df = @(t) that returns RHS derivative for 
%                    a given t
% df_strip_height  - analyticity strip height of df(map_axis2int_direct(t))
% method           - the operator function approximation method. 
%   Available values:
%     'hyperbola'     : method implemented in fcp_op_func_sect_hyp
%     'half-hyperbola': method implemented in fcp_op_func_sect_half_hyp
%     'auto'          : automatic choice based on the maximal value of t
% smoothness_f0    - smoothness parameter of f0 in terms of the operator
%                  powers, i. e. ||A^smoothness_f0 * f0 || < infinity. 
%                  This parameter is optional. (default value: 1)
% smoothness_df    - smoothness parameter of df in terms of the operator
%                  powers, i. e. ||A^smoothness_f0 * f'(t) || < infinity. 
%                  This parameter is optional. (default value: 1)
% 
%
% Parameters u1_res and vu1 are only used when the order of fractional
% derivative in time fcp.frac_order > 1.
%
%OUTPUT:
% fcp_prop_t    - approximated fractional Cauchy propagator evaluated at 
%                 the grid (t,x)
% Sf0           - propagator S operator function structure containing evaluated
%                 resolvents applied to f0
% Sf0_matr_func - matrix of scalar functions associated with the quadrature 
%                 points used to approximate propagator Sf0. 
%                 See cauchy_int_matr_func for details.
% Sdf           - propagator S operator function structure containing evaluated
%                 resolvents applied to df
% Sdf_matr_func - matrix of scalar functions associated with the quadrature 
%                 points used to approximate propagator Sdf. 
%                 See cauchy_int_matr_func for details.
%
% Copyright (C) 2022, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin == 7)
  method = 'auto';
  smoothness_f0 = 1;
  smoothness_df = 1;
elseif (nargin == 8)
  smoothness_f0 = 1;
  smoothness_df = 1;
elseif (nargin == 9)
  smoothness_df = 1;
elseif (nargin == 10)
  % Do nothing
else
  error('Wrong number of input arguments.');
end

% Check and normilize sizes of initial values
if (~isequal(size(f0),size(df(0))))
    error('Two provided f0 and df(0) have incompatible sizes.');
end

if strcmpi(method,'auto')
  if (max(t(:)) > 1)
    method = 'half-hyperbola';
  else
    method = 'hyperbola';
  end
end

%% Function body
% Set whether to use persistent variables or not. If set to true it enables
% the partial resolvent reuse capabilities, so the evaluation of operator 
% function for new values of t can be performed without resolvent 
% re-evaluations, provided that all other input parameters remain the same.
USE_PERSISTENT = false;

% Selecting the operator function implementation
if strcmpi(method,'hyperbola')
  OP_FUNC = @fcp_op_func_sect_hyp;
elseif strcmpi(method,'half-hyperbola')
  OP_FUNC = @fcp_op_func_sect_half_hyp;
else
  error('The user-specified "method" is unknown');
end

% Order of fractional derivative ( used in contour and function specification )
alpha = fcp.frac_order;

%Currently only manual selection for the spectrum parameters of A is
%possible
if strcmp(fcp.spectral_estimator,'manual')
  sphi = fcp.spectral_angle;
  sa0  = fcp.spectral_vertex;
  else
    error('Only manual estimator is implemented so far.');
end

%% Perform data re-use if applicable 
% Persistent variables needed to reuse the previously calculated resolvents
persistent pmethod pfcp pSf0;
persistent pf0_res pf0;

% Check input values and assign previously saved operator structures
% to perform re-evaluation with previously calculated resolvents for each
% of two propagators
op_func_S = [];      % Needed to trigger initialization of op_func_S
if (~isempty(pfcp)) && isequaln(pfcp,fcp) && strcmpi(pmethod,method)
  if isequal(f0_res,pf0_res) && isequal(f0,pf0)
    op_func_S  = pSf0;
  end
elseif (USE_PERSISTENT)
  pfcp = fcp; pmethod = method;
else
  pfcp = []; pmethod = [];
end
%% Define operator structures

% Specification of operator function
% Critical sector is needed for the calculation of contour parameters
critical_sector.shape  = 'sectorial';
critical_sector.vertex = fcp.critical_vertex;
critical_sector.angle  = fcp.critical_angle;

% Specification of operator, sector and resolvent
if (fcp.correction_order > 1)
  if isempty(res_M)
    [~,res_M] = f0_res(0,f0);
  end
  operator.func = @(vu) res_M*vu;
else
  operator.func = [];
end
lap_spectrum.shape  = 'sectorial';
lap_spectrum.vertex = sa0;
lap_spectrum.angle  = sphi;
operator.spectrum   = lap_spectrum;

lap_res.func = f0_res;
lap_res.correction_point = fcp.correction_point;
lap_res.correction_order = fcp.correction_order;
operator.resolvent = lap_res;

%% Contour specification
hyp_contour.type ='hyperbolic';
% Contour parameters are calculated by the corresponding op_func routine
hyp_contour.a0 = NaN;
hyp_contour.ai = NaN;
hyp_contour.bi = NaN;
hyp_contour.b0 = fcp.contour_vertex;

%% Specify and evaluate the propagator operator functions
if isempty(op_func_S)
  % Propagator for f0
  % Make sure that the next function is vectorizable and return the array of
  % size [numel(t),numel(A)]
  op_func_S.func = @(t,A) exp( t .* A).*(A.^(alpha-1));
  op_func_S.contour = hyp_contour;
  op_func_S.func_decay_order = alpha * smoothness_f0;       % Specify the decay order of the integrand (needed to calculate h)
  op_func_S.time = t;                                       % Column vector of times t
  op_func_S.N = fcp.N;                                      % Quadrature discretization parameter
  op_func_S.func_matr = [];                                 % Optional matrix representation of function
  op_func_S.critical_spectrum = critical_sector;            % Critical sector
else  % Re-use previously calculated structure and adjust the time field
  op_func_S.time = t;
  op_func_S.func_matr = [];
end

% Final time 
T  = fcp.final_time;
N  = fcp.N;

% Do the calculations
%1.  Evaluate J(t)*S*f(0)
N0 = ceil(alpha * smoothness_f0 * N / min(1,alpha));
d  = 0.5 * min(pi,(pi - sphi) / alpha) - pi/4;
fcp_prop_t = RLI_approx(@Sf0_func,alpha,t,0,T,N0,d);

% Next two terms are nonzero only if at df(t) != 0
if any(df(t) ~= 0,'all')
  %2. Evaluate int(J*f'(s),s=0..t)
  % r2 = 0;
  N1 = ceil( df_strip_height / d * alpha * smoothness_df * N);
  % N2 is used within RL_df
  N2 = ceil(df_strip_height / d * alpha  * smoothness_df * N / min(1,alpha));
  r21 = quad_SQ_int(@RL_df,0,t(1),N1,df_strip_height);
  if (numel(t) > 1)
    r2 = zeros([numel(r21),numel(t)]);
    r2(:,1) = r21;
    for i = 2:numel(t)
      r2(:,i) = quad_SQ_int(@RL_df,0,t(i),N1,df_strip_height);
    end
  else
    r2 = r21; r21 =[];
  end

  %3. Evaluate int(S(t-s)*J*f(s),s=0..t)
  % r3 = 0;
  N3 = N;
  % N4,N5 are used within fcp_arg_z_J_df
  N4 = ceil(df_strip_height / d * alpha * smoothness_df * N);
  N5 = ceil(df_strip_height / d * alpha * smoothness_df * N / min(1,alpha));
  operator.resolvent.func = df_res;
  op_func_S.func = @(t,A) ones(size(t)) .* (A.^(alpha-1));
  op_func_S.N = N3;
  op_func_S.func_matr = [];
  op_func_S.time = t;

  r3_arg = struct("time_dependent",true,"op_dependent",true,"func",@fcp_arg_z_J_df);
  r3 = OP_FUNC(alpha,op_func_S,operator,r3_arg);

  fcp_prop_t = fcp_prop_t + r2 + r3;
end

%% Local functions
function RL_df_t = RL_df(t)
  RL_df_t = RLI_approx(df,alpha,t,0,T,N2,df_strip_height);
end

function fcp_propS_f0 = Sf0_func(t)
  op_func_S.time = t;
  [fcp_propS_f0,Sf0,Sf0_matr_func] = OP_FUNC(alpha,op_func_S,operator,f0);
  % Save the structures used in the calculation for potential future reuse 
%   if (USE_PERSISTENT) 
%     pf0_res = f0_res; pf0 = f0;
%     pSf0 = Sf0;
%   else
%     pf0_res = []; pf0 = [];
%     pSf0 = [];
%   end
  % Check if the calculated contour encircles the integrand's singularity
  if (Sf0.residue_substraction)
    % Add the term responsible for the residue of the integrand at z=0
    % which is equal to the identity operator
    fcp_propS_f0 = fcp_propS_f0 + repmat(f0,size(t));
  end
end

function G = fcp_arg_z_J_df(t,z)
  fcp_arg_ez_J_df = @(s) exp(z.*(t-s)) .* RLI_approx(df, alpha, s, 0, T, N5, df_strip_height);
  % Prototype: quad_SQ_int(f, a, b, N, d)
  G1 = quad_SQ_int(fcp_arg_ez_J_df, 0, t(1), N4, df_strip_height);
  if (numel(t) > 1)
    G = zeros([numel(G1),numel(t)]);
    G(:,1) = G1;
    for i = 2:numel(t)
      G(:,i) = quad_SQ_int(fcp_arg_ez_J_df,0,t(i),N4,df_strip_height);
    end
  else
    G = G1; 
  end
end
    
end
