function [coeffs, rho] = pick_the_right_coefs(ys, s, Theta, threshold)
% PICK_THE_RIGHT_COEFS computes those coeffs, which would return in a 
% regenerated signal, from which the correlation to the true signal is 
% closest to `threshold`.

[N, M] = size(ys);
assert(size(Theta,2) == N, "Length of coefficient vectors must match the number of basis functions.")
assert(size(Theta,1) == length(s), "Lengths of basis functions must match the length of signal s.")

rhos = zeros(1,M-1);
for i = 1:M
    rr = corrcoef(regenerate_signal(Theta, ys(:,i)), s);
    if isnan(rr(2))
        break
    else
        rhos(i) = rr(2);
    end
end
d = abs(threshold - rhos);
[~, min_idx] = min(d);
coeffs = ys(:, min_idx+1); 
rho = rhos(min_idx);
end
