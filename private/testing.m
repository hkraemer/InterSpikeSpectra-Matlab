clear, clc

N = 50;
period1 = 13;
period2 = 8;
test_tauRR2 = create_single_basis_function(N, period1);
test_tauRR3 = create_single_basis_function(N, period2);


test_tauRR = 0.7 .* test_tauRR2(3,:);
test_tauRR_ = 0.8 .* test_tauRR3(7,:);

test_tauRR = test_tauRR + test_tauRR_;
test_tauRR = test_tauRR + 0.05.*randn(1,N);

%%
s = test_tauRR;

[N, M] = size(s);
if N < M
    s = s';
end
s = s - min(s);
s = s ./ max(s);

%%
threshold = 0.95;
N = length(s);
% get set of basis functions
Theta = generate_basis_functions(N)';
% make LASSO regression
%ys = lasso(Theta, s, 'Intercept', false, 'Standardize', false, 'NumLambda',100);

%%
ys = lasso(Theta, s,'CV',5);



%%
[y, rho] = inter_spike_spectrum(s,.95);

%%

plot(y)