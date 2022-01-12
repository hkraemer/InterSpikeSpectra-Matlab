function [spectrum, rho] = inter_spike_spectrum(s, varargin)
% INTER_SPIKE_SPECTRUM spike decomposition.
%
% computes the inter spike spectrum from the input
% signal `s` (column or line vector). 
%
%    spectrum = INTER_SPIKE_SPECTRUM(X) computes the inter spike spectrum 
%    from the input signal `X` (column or line vector).
%
%    [spectrum, rho] = INTER_SPIKE_SPECTRUM(...) computes the inter spike 
%    spectrum from the input signal `X` (column or line vector) and the 
%    correlation coefficient `rho` of the retransformed decomposed signal 
%    and `X`. 
%       
%    ... = INTER_SPIKE_SPECTRUM(..., Name, Value) specifies further optional 
%    parameters using one or more Name, Value pair arguments.
%
%    Optional name-value-arguments:
%    'threshold'        - (default 0.95) The agreement of the regenerated 
%                         decomposed signal with the true signal. This depends 
%                         on the LASSO regularization parameter lambda. Lambda 
%                         gets adjusted automatically with respect to 'threshold'.
%   'tol'               - (default 1e-3) Allowed tolerance between 'threshold' and 
%                         'rho'.
%   'max_iter'          - (default 15) Determines after how many tried Lambdas 
%                         the algorithm stops and must be an integer between 2 
%                         and 20.
%    'verbose'          - (default true) If true, warning messages enabled.
%
%    Example:
%    Compute the Inter Spike Spectrum of a monochromatic spike signal of spike 
%    period 5.
%
%    N = 50;                                         % signal length
%    period = 5;                                     % spike  period
%    s = zeros(1,N);       
%    s(3:period:end) = 1;
%    s = s * 0.3*randn(1,50);                        % randomize peak heights
%    [spectrum, rho] = inter_spike_spectrum(s);
%    figure
%    subplot(211)
%    plot(1:N,s)
%    title('Signal')
%    grid on
%    subplot(212)
%    plot(1:N,spectrum)
%    title('ISS of signal s')
%    grid on
% 
%    Further reading:
%    H. K. Kraemer et al., … 2022
%
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
narginchk(1,9)
nargoutchk(1,2)

% Default values
threshold = 0.95;
tol = 1e-3;
max_iter = 15;
verbose = true;

% required and optional arguments
p = inputParser;

validScalarPosNum0 = @(x) isnumeric(x);
validScalarPosNum1 = @(x) isnumeric(x) && isscalar(x) && (x >= 0.8) && (x <= 1);
validScalarPosNum2 = @(x) isnumeric(x) && isscalar(x) && (x >= 1e-5) && (x <= 1);
validScalarPosNum3 = @(x) isnumeric(x) && isscalar(x) && (x <= 20) && (x > 1);
validScalarPosNum4 = @(x) islogical(x);

addRequired(p, 's', validScalarPosNum0);
addParameter(p, 'threshold', threshold, validScalarPosNum1);
addParameter(p, 'tol', tol, validScalarPosNum2);
addParameter(p, 'max_iter', max_iter, validScalarPosNum3);
addParameter(p, 'verbose', verbose, validScalarPosNum4);

% parse input arguments
parse(p,s,varargin{:})

% assign variables with the resulting argument input
s = p.Results.s;
threshold = p.Results.threshold;
tol = p.Results.tol;
max_iter = p.Results.max_iter;
verbose = p.Results.verbose;

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
% make LASSO regression
[spectrum, rho] = compute_spectrum_according_to_threshold(s, Theta, threshold, tol, max_iter, verbose);

end