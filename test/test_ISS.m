% Test InterSpikeSpectra functionality

clear, clc

cd ../private/
% preconditions
rng(1234)
N = 50;
period1 = 13;
period2 = 8;
test_tauRR2 = create_single_basis_function(N, period1);
test_tauRR3 = create_single_basis_function(N, period2);
test_tauRR = 0.7 .* test_tauRR2(3,:);
test_tauRR_ = 0.8 .* test_tauRR3(7,:);
test_tauRR = test_tauRR + test_tauRR_;
test_tauRR = test_tauRR + 0.05.*randn(1,N);

cd ../

%% Test 1: 2 Spikes + measurement noise
[~, rho] = inter_spike_spectrum(test_tauRR, 0.85);
assert(0.83 <= rho && rho < 0.89)

[spectrum, rho] = inter_spike_spectrum(test_tauRR, 0.99);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.1;
peak_idxs = max_idx(t_idx);

assert(0.989 <= rho && rho < 0.991)
assert(length(peak_idxs) == 2)
assert(peak_idxs(1) == period2)
assert(peak_idxs(2) == period1)

%% Test 2: randomized peak heights
rng(1234)
test_tauRR = abs(randn(1,N)) .* test_tauRR2(1,:);
[spectrum, ~] = inter_spike_spectrum(test_tauRR);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.01;
peak_idxs = max_idx(t_idx);

assert(sum(rem(peak_idxs, period1)) == 0)