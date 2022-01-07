function [spectrum, rho] = inter_spike_spectrum(varargin)
% INTER_SPIKE_SPECTRUM computes the inter spike spectrum from the input
% signal `s` (column or line vector). 
%
%   spectrum [, corr] = inter_spike_spectrum(s [,correlation_threshold, tol, maxlambdas])
%
% The optional second input argument `correlation_threshold` determines how
% well the retransformed decomposed signal correlates with the given input 
% `s` (linear Pearson correlation coefficient, Default is 
% `correlation_threshold=0.95`). The regularization parameter 
% Lambda is then automatically adjusted by tracking the matching within a 
% given tolerance `tol`, which is the third optional input (Default is 
% `tol=1e-3`). The fourth optional input `maxlambdas` determines after how
% many tried Lambdas the algorithm stopps (Default is 15).

% The second optional output `corr` is the obtained correlation coefficient 
% of the retransformed decomposed signal and `s` and should match 
% `correlation_threshold`.
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
narginchk(1,4)
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
assert(threshold >= 0.8, "Optional input `correlation_threshold` must be a value in the interval [0.8, 1]")
assert(threshold <=1, "Optional input `correlation_threshold` must be a value in the interval [0.8, 1]")

try
    tol = varargin{3};
catch
    tol = 1e-3;
end
assert(isnumeric(tol), "Optional input `tol` must be numeric.")
assert(tol > 1e-5, "Optional input `tol` must be a value in the interval (1e-5, 1]")
assert(tol <=1, "Optional input `tol` must be a value in the interval (1e-5, 1]")


try
    maxlambdas = varargin{4};
catch
    maxlambdas = 15;
end
assert(rem(maxlambdas,1)==0, "Optional input `maxlambdas` must be an integer in the interval (0, 100].")
assert(maxlambdas <= 100, "Optional input `maxlambdas` must be an integer in the interval (0, 100].")
assert(maxlambdas > 0, "Optional input `maxlambdas` must be an integer in the interval (0, 100].")

% normalize time series
s = (s - mean(s)) ./ std(s);
s = s - min(s);
s = s ./ max(s);

%% Compute spectrum

N = length(s);
% get set of basis functions
Theta = generate_basis_functions(N)';
% make LASSO regression
[spectrum, rho] = compute_spectrum_according_to_threshold(s, Theta, threshold, tol, maxlambdas);

end