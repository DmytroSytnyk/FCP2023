function [fcp_prop_t,S,S_matr_func,P,P_matr_func] = fcp_cauchy_hom_appr(t,fcp,u0_res,u1_res,u0,u1,method)

% FCP_CAUCHY_HOM_APPR approximate solution of homogeneous fractional Cauchy problem.
%
%INPUT:
% t       - the array time-points for the solution approximation
% fcp     - fractional Cauchy problem (FCP) structure that defines problem 
%           and operator parameters. Structure fcp is assumed to have the 
%           following fields:
%    frac_order      : fractional derivative order parameter
%    spectral_angle  : actual spectral angle
%    spectral_vertex : actual spectral vertex
%    critical_angle  : critical angle
%    critical_vertex : critical vertex
%    contour_vertex  : the prescribed value of contour at zero
% u0_res  - function implementing the evaluation of resolvent applied to vu0
% u1_res  - function implementing the evaluation of resolvent applied to vu1
%   For FCP the resolvent has to be defined as
%                   R(z) = z^(alpha-1)*(z^alpha*I + A)^(-1)
%   whence the minus sign in the definition of the contour below.
%   So in order to reuse the resolvent functions that implement convenient
%   resolvent implementation  R(z) = (z*I - A)^(-1) we will define
%   adhoc function wrapper. This will allow to freely reuse the resolvent
%   implementations across different problems.
% u0     - first initial value
% u1     - second initial value
%   u0,u1  might be a vector or structure with two fields:
%   time_dependent: Flag indicating whether the func is time-dependent
%   func          : Func having prototype @(t) (time_dependent flag = true)
%                   or @() (time_dependent = false)
% method - the operator function approximation method. Available values:
%   'hyperbola'      : method implemented in fcp_op_func_sect_hyp
%   'half-hyperbola' : method implemented in fcp_op_func_sect_half_hyp
%   'auto'           : automatic choice based on the maximal value of t
%
% Parameters u1_res and vu1 are only used when the order of fractional
% derivatime in time fcp.frac_order > 1.
%
%OUTPUT:
% fcp_prop_t  - approximated fractional Cauchy propagator evaluated at the
%               grid (t,x)
% S           - progator S operator function structure containing evaluated
%               resolvents
% S_matr_func - matrix of scalar functions associated with the quadrature
%               points used to approximate propagator S.
%               See cauchy_int_matr_func for details.
% P           - progator P operator function structure containing evaluated
%               resolvents
% P_matr_func - matrix of scalar functions associated with the quadrature
%               points used to approximate propagator P.
%               See cauchy_int_matr_func for details.
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input parameters
if (nargin == 5)
  if (fcp.frac_order > 1)
    error('Two initial values must be supplied if the fractional order is greater than one.');
  end
  method = 'auto';
elseif (nargin == 6)
  method = 'auto';
elseif (nargin == 7)
  % Do nothing
else
  error('Wrong number of input parameters');
end

% Parse initial values and do sanity checks
u0 = parse_initial_value(u0);
if (u0.time_dependent)
  vu0 = u0.func(t(1));
else
  vu0 = u0.func();
end
u1 = parse_initial_value(u1);
if (u1.time_dependent)
  vu1 = u1.func(t(1));
else
  vu1 = u1.func();
end

% Check and normalize sizes of initial value vectors
if (size(vu0) ~= size(vu1))
  if (numel(vu1) == 1)
    vu1 = vu1(1)*ones(size(vu0));
  elseif (numel(vu0) == 1)
    vu0 = vu0(1)*ones(size(vu1));
  else
    error('Two provided initial values should have compatible sizes');
  end
end

if strcmpi(method,'auto')
  if (max(t(:))>1)
    method = 'half-hyperbola';
  else
    method = 'hyperbola';
  end
end

%% Function body
% Set whether to use persistent variables or not
USE_PERSISTENT = true;
smoothness_u0 = 1;
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

% We use a black box approach to estimate the spectrum parameters of A:
% 1. Calculate the matrix discretization of differential operator
% 2. Find the extreme eigenvalues of that matrix
if strcmp(fcp.spectral_estimator, 'manual')
  sphi = fcp.spectral_angle;
  sa0 = fcp.spectral_vertex;
  % else
  %   error('Only manual estimator is implemented so far.');
elseif strcmp(fcp.spectral_estimator, 'eigs_real')
  [res_approx,res_M] = u0_res(0,vu0);
  %    res_M = full(res_M);
  opts.tol = 1e-3;
  opts.maxit = 1e4;
  [~,lambda,flag] = eigs(res_M, 1,'sr', opts);
  if flag ~=0 || imag(lambda)>0
    error('Failed to estimate the spectrum lowest real eigenvalue. Try eigs_complex.');
  end
  sa0 = real(lambda);
  sphi = 0;
elseif strcmp(fcp.spectral_estimator, 'eigs_complex')
  [res_approx,res_M] = u0_res(0,vu0);
  %    res_M = full(res_M);
  opts.tol = 1e-5;
  opts.maxit = 1e4;
  [~,lambda,flag] = eigs(res_M, 1,'sr', opts);
  if flag ~=0
    error('Failed to estimate the lowest real eigenvalue');
  end
  lambda_r = real(lambda);
  [~,lambda,flag] = eigs(res_M, 1,'si', opts);
  if flag ~=0
    error('Failed to estimate the lowest imaginary eigenvalue');
  end
  lambda_si = imag(lambda);
  [~,lambda,flag] = eigs(res_M, 1,'li', opts);
  if flag ~=0
    error('Failed to estimate the largest imaginary eigenvalue');
  end
  lambda_li = imag(lambda);
  lambda_i = max(lambda_si, lambda_li);
  sphi = pi/4;
  sa0 = lambda_r - lambda_i*tan(sphi);
elseif strcmp(fcp.spectral_estimator, 'eig')
  %   error('Not implemented');
  [~,res_M] = u0_res(0,vu0);
  lambda = eig(full(res_M));
  % Make sure that oe object defines sector
  fcp.spectral_d = 0;
  [~,~,~,~,d,d0,b_s,d_s,sphi] = cont_trunc_sect_est_pars(lambda, fcp.spectral_vertex,fcp.spectral_d,fcp.spectral_angle,fcp.spectral_angle,fcp.critical_angle);
  if (d_s > 0)
    error('Internal error in sectorial spectrum estimator');
  else
    sa0 = b_s;
  end
end

%% Perform data reuse if applicable
% Persistent variables needed to reuse the previously calculated resolvents
persistent pmethod pfcp pS pP;
persistent pu0_res pu1_res pvu0 pvu1;

% Check input values and assign previously saved operator structures
% to perform re-evaluation with previously calculated resolvents for each
% of two propagators
op_func_S = [];      % Needed to trigger initialization of op_func_S
op_func_P = [];      % Needed to trigger initialization of op_func_P
if (~isempty(pfcp)) && isequaln(pfcp,fcp) && strcmpi(pmethod,method)
  if isequal(u0_res,pu0_res) && isequal(vu0,pvu0)
    op_func_S  = pS;
  end
  if (alpha > 1) && isequal(u1_res,pu1_res) && isequal(vu1,pvu1)
    op_func_P  = pP;
  end
elseif (USE_PERSISTENT)
  pfcp = fcp; pmethod = method;
else
  pfcp = []; pmethod = [];
end

%% Define operator data structures
% Specification of operator function
% Critical sector is needed for the calculation of contour parameters
critical_sector.shape = 'sectorial';
critical_sector.vertex = fcp.critical_vertex;
critical_sector.angle = fcp.critical_angle;

% Specification of operator, sector and resolvent
if (fcp.correction_order > 1)
  if isempty(res_M)
    [~,res_M] = u0_res(0,vu0);
  end
  operator.func = @(vu) res_M*vu;
else
  operator.func = [];
end
lap_spectrum.shape = 'sectorial';
lap_spectrum.vertex = sa0;
lap_spectrum.angle = sphi;
operator.spectrum = lap_spectrum;

lap_res.func = u0_res;
lap_res.correction_point = fcp.correction_point;
lap_res.correction_order = fcp.correction_order;
operator.resolvent = lap_res;

%% Contour specification
hyp_contour.type ='hyperbolic';

% The contour parameters are calculated by the corresponding fcp_op_func routine
hyp_contour.a0 = NaN;
hyp_contour.ai = NaN;
hyp_contour.bi = NaN;
hyp_contour.b0 = fcp.contour_vertex;

if isempty(op_func_S)
  %% Specify and evaluate the propagator operator functions
  % Propagator for alpha <=1
%   op_func_S = struct([]);
  op_func_S.contour = hyp_contour;
  
  % Make sure that the next function is vectorizable
  op_func_S.func = @(t,A) exp( t .* A).*(A.^(alpha-1));
  op_func_S.func_decay_order = alpha * smoothness_u0;       % Specify the decay order of the integrand (needed to calculate h)
  op_func_S.time = t;                                       % Column vector of times t
  op_func_S.N = fcp.N;                                      % Quadrature discretization parameter
  op_func_S.func_matr = [];                                 % Optional matrix representation of function
  op_func_S.critical_spectrum = critical_sector;            % Critical sector
else  % Re-use previously calculated structure and adjust the time field
  op_func_S.time = t;
  op_func_S.func_matr = [];
end
% Do the calculations
[fcp_prop_t,S,S_matr_func] = OP_FUNC(alpha,op_func_S, operator, u0);
  % Save the structures used in the calculation for potential future reuse 
if (USE_PERSISTENT) && (~isempty(S_matr_func)) 
  pu0_res = u0_res; pvu0 = vu0;
  pS = S;
else
  pu0_res = []; pvu0 = [];
  pS = [];
end
% Check if the calculated contour encircles the signularity of integrand
if (S.residue_substraction)
  % Add the term responsible for the residue of the integrand at z=0
  % which is equal to the identity operator
  if (u0.time_dependent)
    fcp_prop_t = fcp_prop_t + horzcat(vu0,u0.func(t(2:end)));
    %NOTE! u0.func should return the proper size array for vector input
  else
    fcp_prop_t = fcp_prop_t + repmat(vu0,size(t));
  end
end
if (alpha > 1)
  % Additional propagator (int(S(s),s=0..t)) needed only if alpha > 1
  % Change the operator resolvent to match the initial values of this
  % propagator
  lap_res.func = u1_res;
  lap_res.correction_point = fcp.correction_point;
%   lap_res.correction_order = fcp.correction_order;
  lap_res.correction_order = 0;
  operator.resolvent = lap_res;
  
  if isempty(op_func_P)
%     op_func_P = struct([]);
    % Reuse the contour parameters
    op_func_P.contour = op_func_S.contour;
    % Make sure that the next function is vectorizable
    op_func_P.func = @(t,A) (exp( t .* A)).*(A.^(alpha-2));
    op_func_P.func_decay_order = 1;                         % Specify the decay order of the integrand (needed to calculate h)
    op_func_P.time = t;                                     % Column vector of times t
    op_func_P.N = fcp.N;                                    % Quadrature discretization parameter
    op_func_P.func_matr = [];                               % Optional matrix representation of function
    op_func_P.critical_spectrum = critical_sector;          % Critical sector
  else  % Re-use previously calculated structure and adjust the time field
    op_func_P.time = t;
  end
  % Do the calculations
  [fcp_prop_Pt,P,P_matr_func] = OP_FUNC (alpha,op_func_P, operator, vu1);
  % Save the structures used in the calculation for potential future reuse 
  if (USE_PERSISTENT)
    pu1_res = u1_res; pvu1 = vu1;
    pP = P;
  else
    pu1_res = []; pvu1 = [];
    pP = [];
  end
  % Check if the calculated contour encircle residue of the integrand
%   if (P.residue_substraction)
%     % Add the term responsible for the residue of the integrand at z=0
%     % which is equal to the identity operator
%     if u1.time_dependent
%       fcp_prop_Pt = fcp_prop_Pt + repmat(vu1,size(t));
%     else
%       fcp_prop_Pt = fcp_prop_Pt + horzcat(vu1,u1.func(t(2:end)));
%     end
%   end
  fcp_prop_t = fcp_prop_t + fcp_prop_Pt;
end

%% Nested functions
  function arg_out = parse_initial_value(arg)
    arg_out = struct("time_dependent",false,"op_dependent",false, "func", []);
    if isa(arg,"double")
      arg_out.func = @() arg;
    elseif isa(arg,"struct")
      if isfield(arg, "time_dependent")
        arg_out.time_dependent =  arg.time_dependent;
      end
      if isfield(arg, "op_dependent")
        arg_out.op_dependent =  arg.op_dependent;
      end
      if isfield(arg, "func")
        arg_out.func =  arg.func;
      end
    else
      error('Argument has wrong type. Consider using a matrix or the struct("time_dependent", "op_dependent", "func").')
    end
  end
end

