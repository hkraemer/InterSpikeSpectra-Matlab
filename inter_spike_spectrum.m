function [spectrum, rho] = inter_spike_spectrum(s, varargin)
% INTER_SPIKE_SPECTRUM spike decomposition.
%
%    S=INTER_SPIKE_SPECTRUM(X) computes the inter spike spectrum S
%    from the input signal X (column or line vector). 
%
%    [S, R]=INTER_SPIKE_SPECTRUM(...) computes the inter spike 
%    spectrum S from the input signal X (column or line vector) and  
%    the correlation coefficient R between the retransformed 
%    decomposed signal and X. 
%       
%    ...=INTER_SPIKE_SPECTRUM(..., Name, Value) specifies further optional 
%    parameters using one or more (Name, Value) pair arguments.
%
%    Optional name-value-arguments:
%    'method'    - (default "lasso") The method for sparse regression. Pick 
%                  either "lasso" or "STLS" (sequential thresholded least squares)
%    'logit'     - (default true) Whether or not a logistic regression
%                  should be made or not.
%    'threshold' - (default 0.99) The agreement of the regenerated 
%                  decomposed signal with the true signal. This depends 
%                  on the LASSO regularization parameter lambda. Lambda 
%                  gets adjusted automatically with respect to 'threshold'.
%    'Alpha'     - (default is 1, which corresponds to LASSO) The wheight
%                  between LASSO and Ridge regression when 'method' == "lasso". 
%                  1 corresponds to LASSO and 0 to pure Ridge regression. Values 
%                  between 0 and 1 correspond to an elastic net.
%    'tol'       - (default 1e-3) Allowed tolerance between 'threshold' and R.
%    'max_iter'  - (default 15) Determines after how many tried lambdas the
%                  algorithm stops and must be an integer between 2 and 20.
%    'verbose'   - (default true) If true, warning messages enabled.
%
%    Example:
%    % Compute the Inter Spike Spectrum of a monochromatic spike signal of spike 
%    period 5.
%
%    N = 100;                                         % signal length
%    period = 5;                                     % spike  period
%    s = zeros(1,N);       
%    s(3:period:end) = 1;
%    s = s + 0.05.*randn(1,N);                        % randomize peak heights
%    spectrum = inter_spike_spectrum(s);
%    subplot(211)
%    plot(s)
%    title('Signal')
%    grid on
%    subplot(212)
%    area(spectrum)
%    title('ISS of signal s')
%    grid on
% 
%    Further reading:
%    H. K. Kraemer et al., … 2022
%

% Copyright (c) 2022-
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Physics, University of Potsdam, Germany
% http://www.physik.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% Check input & output
narginchk(1,15)
nargoutchk(1,2)

% Default values
threshold = 0.99;
tol = 1e-3;
max_iter = 15;
verbose = true;
method = "lasso";
logit = true;
Alpha = 1;

% required and optional arguments
p = inputParser;

validScalarPosNum0 = @(x) isnumeric(x);
validScalarPosNum1 = @(x) isnumeric(x) && isscalar(x) && (x >= 0.8) && (x <= 1);
validScalarPosNum2 = @(x) isnumeric(x) && isscalar(x) && (x >= 1e-5) && (x <= 1);
validScalarPosNum3 = @(x) isnumeric(x) && isscalar(x) && (x <= 20) && (x > 1);
validScalarPosNum4 = @(x) islogical(x);
validScalarPosNum5 = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1);
validString = @(x) isstring(x) && (strcmp(x,"lasso") || strcmp(x,"STLS"));

addRequired(p, 's', validScalarPosNum0);
addParameter(p, 'threshold', threshold, validScalarPosNum1);
addParameter(p, 'method', method, validString);
addParameter(p, 'tol', tol, validScalarPosNum2);
addParameter(p, 'max_iter', max_iter, validScalarPosNum3);
addParameter(p, 'verbose', verbose, validScalarPosNum4);
addParameter(p, 'logit', logit, validScalarPosNum4);
addParameter(p, 'Alpha', Alpha, validScalarPosNum5);

% parse input arguments
parse(p,s,varargin{:})

% assign variables with the resulting argument input
s = p.Results.s;
threshold = p.Results.threshold;
tol = p.Results.tol;
max_iter = p.Results.max_iter;
verbose = p.Results.verbose;
method = p.Results.method;
logit = p.Results.logit;
Alpha = p.Results.Alpha;

[N, M] = size(s);
if N < M
    s = s';
end
assert(M == 1 || N == 1, "Input signal must be a column or row vector.")

% normalize time series
s = (s - mean(s)) ./ std(s);
s = s - min(s);
s = s ./ max(s); 

%% Compute spectrum
N = length(s);
% get set of basis functions
Theta = generate_basis_functions(N)';
% make sparse regression
if strcmp(method,"lasso")
    if logit
        [spectrum, rho] = compute_spectrum_according_to_threshold_lasso_logit(s, Theta, threshold, tol, max_iter, Alpha, verbose);
    else
        [spectrum, rho] = compute_spectrum_according_to_threshold_lasso(s, Theta, threshold, tol, max_iter, Alpha, verbose);
    end
elseif strcmp(method,"STLS")
    [spectrum, rho] = compute_spectrum_according_to_threshold_stls(s, Theta, threshold, tol, max_iter, verbose);
end
