function [op_func_appr, op_func_mod, matr_func] = fcp_op_func_sect_hyp(alpha, op_func, operator, arg)

%FCP_OP_FUNC_SECT_HYP - evaluates the action of sectorial operator function using hyperbolic contour
%
% This function numerically evaluates the action of
% function arising from the representation of progators to fractional
% Cauchy problem (analytic function bounded in the
% left part of the complex plane) of operator A on the vector vu
% at a given time t (t can be vector).
%
%INPUT:
% alpha      - fractional order parameter
% op_func    - structure that defines operator function
% operator   - structure that defines operator
% arg        - argument that operator function acts upon. It might be
%              vector or structure with 3 fields:
%  time_dependent: Flag indicating whether the func is time-dependent 
%  op_dependent  : Flag indicating whether the func is time-dependent 
%  func          : Func having prototype:
%                 @()    (if time_dependent = false, op_dependent = false)
%                 @(t)   (if time_dependent = true, op_dependent = false)
%                 @(t,z) (if time_dependent = true, op_dependent = true)
%                 z dependent argumens are not supported
%
%OUTPUT:
% op_func_appr - the result of the approximate operator function evaluation
% op_func_mod  - modified opertor function structure
% matr_func    - the seqence of values of scalar version of the input 
%                function evaluate at the quadrature points. 
%                This is an optional output parameter. Needed for the
%                purpose of computation reuse
%
%EXAMPLE:
% op_func_appr = int(exp(z*t)*z^(alpha-1)*(z^alpha*I + A)^(-1)vu, z \in Gamma)
%
% Integration contour is the oposite of hyperbolic contour Gamma: -z where
% z = a0 + ai*cosh(x) - I * bi*sinh(x)
% The operator A is assumed to be sectorial
%   (spectrum is contained within a sector)
% sa0, sphi -  vertex and angle of spectral sector
% ca0 - real value of critical line (ca0 < sa0)
% N - time discretization parameter
% vu - the given vector
% f_res - function handle for resolvent calculation:
% f_res = @(z,vu)
% a0 = ca0;
% ai = (sa0^(1/alpha)-a0)*(cos(d)-sin(d)*tan(phis/alpha));
% bi = (sa0^(1/alpha)-a0)*(cos(d)*tan(phis/alpha)+sin(d));
%
% The parameters are expected to have the following structure
%%% Operator function specification
%op_funct.func = @(t,A) exp( - t .* A); % Make sure it is vectorizable
%op_func.time = [0;2;1];  % Column vector of times t
%op_func.N = 1024;  % Quadrature discretization parameter
%op_func.func_matr = []; %Optional matrix representation of function
%% If the parameter func_matr is provided and compatible with N it
%% overrides parameters func and time
%
%%Type of spectral estimator. Possible values are:
%% 'eig'          - using MATLAB eig function (exact eigenvalues). Use
%%                  spectral_phi ONLY then spectral is consist of one
%%                  eigenvalue.
%% 'eigs_real'    - using MATLAB eig function (exact eigenvalues)
%% 'eigs_complex' - using MATLAB eig function (exact eigenvalues)
%% 'manual'       - user provided parameters (spectral_vertex, spectral_phi)
%
%
%op_func.spectral_estimator = 'eig';
%op_func.critical_spectrum.shape = 'sectorial';
%op_func.critical_spectrum.angle = pi/2;
%op_func.critical_spectrum.vertex = NaN;
%
%%% Operator specification
% Function to calculate action of A on vu
%operator.func = laplace(vu)
%lap_spectrum.shape = 'sectorial';
%lap_spectrum.vertex = -1;
%lap_spectrum.angle = pi/4;
%operator.spectrum = lap_spectrum;
%lap_res.func = laplace_resolvent(z, vu);
%lap_res.correction_point = NaN;
%lap_res.correction_order = 1;
%operator.resolvent = lap_res;
%
%%% Contour specification
%hyp_contour.type ='hyperbolic';
%
%% The contour parameters are calculated by the sec_op_func routine
%hyp_contour.a0 = NaN;
%hyp_contour.ai = NaN;
%hyp_contour.bi = NaN;
%op_func.contour = hyp_contour;
%
%% The argument must be a column vector
%vu = eye(10,1)
%
% Copyright (C) 2022-2025, Dmytro Sytnyk. All rights reserved. 

%% Parse input arguments
if isa(arg,'double')
  vu = arg;
  arg = struct("time_dependent",false,"op_dependent",false, "func", @()vu);
elseif isa(arg,'struct')
  % Do nothing
else
  error('Argument has wrong type. Consider use array of doubles or structure("time_dependent","op_dependent","func").')
end

if ~strcmpi(op_func.critical_spectrum.shape,"sectorial")
  error('Spectral estimation strategy for such critical spectrum shape is not implemented');
end

%% Function body
op_func_mod = op_func;

sphi = operator.spectrum.angle;
% The propagator has singularity at zero which should be treated the same
% way as a part of spectrum hence we override operator.spectrum.vertex with
% zero
% sa0  = min([sign(operator.spectrum.vertex)*abs(operator.spectrum.vertex)^(1/alpha),0]);
sa0 = 0;
% The critical parameters are specified by the behaviour of the scalar part
% of the operator function, hence no alterations are done here
cphi = op_func.critical_spectrum.angle;
ca0  = op_func.critical_spectrum.vertex;
if (ca0 > 0) 
  error(['The integration contour for the propagator should encircle zero.\n',
         'Make sure that the critical vertex < 0.']);
end
% Calculate the parameters of hyperbolic contour
[a0,ai,bi,d] = cont_hyperbolic_get_pars_frac(alpha,sa0,sphi,ca0,cphi);
% For this function implementation resolvent has to be defined as
% as R(z) = z^(alpha-beta)*(z^alpha*I + A)^(-1)
% with beta = 1 for the first propagator (alpha <= 1 ) and beta =2 for the
% second propagator (alpha > 1)
% whence the alpha power and the in the argument of resolvent evaluation 
% function below
resolvent_func = @(z,arg) -operator.resolvent.func(-z,arg);
% and the minus sign in the definition of the contour below
scalar_part_cont_func = @(z) - cont_hyperbolic(z, a0, ai, bi);             % Used for the evaluation of scarap part
% Note that the resolvent part's contour function is also defined and in
% principle might differ from scalar_part_cont_func

% Update the output structure
op_func_mod.contour.a0 = a0;
op_func_mod.contour.ai = ai;
op_func_mod.contour.bi = bi;
% op_func_mod.contour.func = contour_func;

% Determine the type of the argument and, if function, sample one argument 
% value for later use 
if (arg.time_dependent)
  if (arg.op_dependent)
    vu = arg.func(op_func.time(1),scalar_part_cont_func(0));
  else
    vu = arg.func(op_func.time(1));
  end
else
  if (arg.op_dependent)
    vu = arg.func(scalar_part_cont_func(0));
  else
    vu = arg.func();
  end
end

% Calculation of quadrature parameters
N = op_func.N;
if isfield(op_func,'func_decay_order')
  quad_h = quad_SQ_get_h(N, d, op_func.func_decay_order);
else
  assert(false,'Function decay order was not provided. Please supply "func_decay_order" field.')
  quad_h = quad_SQ_get_h(N, d, alpha);
end


% Calculation of the resolvent correction point
cor_order = operator.resolvent.correction_order;
cor_point = operator.resolvent.correction_point;
op_func_mod.residue_substraction = false;
if (cor_order > 0) 
%   if (~arg.op_dependent)
    if (cor_order > 1)
      error('Correction order higher than one is not implemented yet');
    end
    if (~isnan(cor_point))
      %   fprintf('Enforcing the correction point: %e', cor_point);
      if (real(cor_point - scalar_part_cont_func(-1i*d)) < 0) ...
          && (cor_order > 0) && (~arg.op_dependent)
        %     warning(['The correction point lies inside the integration contour.', ' Additional residue term will be added to approximation.']);
        op_func_mod.residue_substraction = true;
      end
    else
      % WARNING! There is no guarantee that the corrected resolvent will
      %          be small enough with the following choice of cor_point
      % TODO: This may require tweaking
      cor_point = scalar_part_cont_func(-1i*d) + 1;
    end
%   else
%     error(['Apply correction to operator-dependent arguments are ill-defined.', ...
%       'Please do so manually or make sure that the arg has op_dependent field set to false.']);
%   end
end
nvu = numel(vu);

z_idx = (0:N);
nz = numel(z_idx);

res_array = [];
if (~arg.time_dependent) && (~arg.op_dependent) ...
    && isfield(op_func,'res_array') && (~isempty(op_func.res_array)) 
  fprintf('The array of calculated resolvent values was provided. Attempting to reuse...');
  if (size(op_func.res_array) == [nvu nz])
    %     freevaluate = true;
    res_array = op_func.res_array;
    fprintf('sucess.\n');
  else
    fprintf('failed.\n The provided array`s size is incompatible with other function parameters.\n');
  end
end

% At this point calculations can be parallelized
% Check if the provided matrix representation is compatible with z range
nt = numel(op_func.time);
if (size(op_func.func_matr) == [nt, nz])
  matr_func = op_func.func_matr;
else
  t = reshape(op_func.time, nt, 1);
  % The next step can be quite time consuming. We can overlap it with
  % evaluation of resolvents.
  matr_func = op_func.func(t, scalar_part_cont_func(quad_h * z_idx));
  if ~isequal(size(matr_func),[nt, nz])
    error(sprintf(['The output of the operator function has size incompatible to size of its input arguments.\n',...
           'Make sure that the provided function is able to handle arrays and its output size\n',...
           'depends on the sizes of both inputs.']));
  end
end

if (arg.time_dependent)
  % For time dependent argument we have to perform calculation for each t_i
  % separately because resolvents can not be reused for different times
  op_func_appr = zeros(numel(vu),nt);
%   if (arg.op_dependent)
%   else
%   end
  if (z_idx(1) == 0) 
    % Half-contour evaluation is suitable for real valued resolvents only
    % To properly use half of the contour we need to divide f(z(0)) by 2
    % and, after the evaluation is done take twice real part of the result
    % as an output
    matr_func(:,1) = 0.5*matr_func(:,1);
  end  
  % parfor i = 1:numel(t)
  for i = 1:numel(t)
    if (arg.op_dependent)
      % We need to iterate trough indices and sum up manually since the
      % argument can not be evaluated fo multiple z_k at once 
      r = zeros(size(vu));
      for k = 1:numel(z_idx)
        vu = arg.func(t(i), scalar_part_cont_func(z_idx(k)*quad_h));
        r = r + cauchy_int_matr_func(operator.func, resolvent_func, @resolvent_part_cont_func, matr_func(i,k), vu, cor_order, cor_point, quad_h, res_array, z_idx(k));
      end
      op_func_appr(:,i) = r;
    else
      vu = arg.func(t(i));
      %NOTE! At present the code does not support reevaluation for time
      % dependent arguments, hence only one output parameter is required
      op_func_appr(:,i) = cauchy_int_matr_func(operator.func, resolvent_func, @resolvent_part_cont_func, matr_func(i,:), vu, cor_order, cor_point, quad_h, res_array, z_idx);
    end
  end
  if (z_idx(1) == 0)
    op_func_appr = 2*real(op_func_appr);
  end
  %   op_func_mod.res_array = res_array_out;
  % At present we disable res_array return for time dependent arguments
  op_func_mod.res_array = [];
else
  % Argument is not time dependent, hence we can evaluate the operator
  % function for all t's in one call to cauchy_int_matr_func
  if (z_idx(1) == 0) 
    % Half-contour evaluation for suitable for real valued resolvents only
    % To properly use half of the contour we need to divide f(z(0)) by 2
    matr_func(:,1) = 0.5*matr_func(:,1);
    [op_func_appr,res_array_out] = cauchy_int_matr_func(operator.func, resolvent_func, @resolvent_part_cont_func, matr_func, vu, cor_order, cor_point, quad_h, res_array, z_idx);
    op_func_appr = 2*real(op_func_appr);
  else  
    % Full-contour evaluation suitable for general complex resolvents  
    [op_func_appr,res_array_out] = cauchy_int_matr_func(operator.func, resolvent_func, @resolvent_part_cont_func, matr_func, vu, cor_order, cor_point, quad_h, res_array, z_idx);
  end
  op_func_mod.res_array = res_array_out;
end

%% Nested functions 
  function [cont,dcont] = resolvent_part_cont_func(x)
    [cont,dcont]  = cont_hyperbolic(x, a0, ai, bi);
    cont =  (-cont).^alpha;      % Used for evaluation of resolvent and correction
    dcont = -dcont;              % Derivative should be unaffected by alpha
  end

end


