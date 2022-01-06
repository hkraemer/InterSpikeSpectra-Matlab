function [spectrum, rho] = inter_spike_spectrum(varargin)
% INTER_SPIKE_SPECTRUM computes the inter spike spectrum from the input
% signal `s` (column or line vector). 
%
%   spectrum [, corr] = inter_spike_spectrum(s [,correlation_threshold])
%
% The optional second input argument `correlation_threshold` determines how
% well the retransformed decomposed signal correlates with the given input 
% `s` (linear Pearson correlation coefficient). The regularization parameter 
% Lambda is then automatically adjusted. Default is 
% `correlation_threshold=0.95`. The second optional output `corr` is the 
% obtained correlation coefficient of the retransformed decomposed signal 
% and `s` and should match `correlation_threshold`.
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
narginchk(1,2)
nargoutchk(1,2)

s = varargin{1};
[N, M] = size(s);
if N < M
    s = s';
end
assert(M == 1 || N == 1, "Input signal must be a column or row vector.")

try
    threshold = varargin{2};
catch
    threshold = 0.95;
end
assert(isnumeric(threshold), "Optional input `correlation_threshold` must be numeric.")
assert(threshold > 0, "Optional input `correlation_threshold` must be a value in the interval (0, 1]")
assert(threshold <=1, "Optional input `correlation_threshold` must be a value in the interval (0, 1]")

% normalize time series
s = s - min(s);
s = s ./ max(s);

%% Compute spectrum

N = length(s);
% get set of basis functions
Theta = generate_basis_functions(N)';
% make LASSO regression
%ys = lasso(Theta, s, 'Intercept', false, 'Standardize', false, 'NumLambda',100);
ys = lasso(Theta, s, 'NumLambda', 100);

[y, rho] = pick_the_right_coefs(ys, s, Theta, threshold);

% pool the same frequencies
spectrum = zeros(1,N);
cnt = 1;
for i = 1:N
    occs = sum(y(cnt:cnt+i-1)>0);
    if occs>0
        spectrum(i) = sum(y(cnt:cnt+i-1))/occs;
    else
        spectrum(i) = 0;
    end
    cnt = cnt + i;
end
end