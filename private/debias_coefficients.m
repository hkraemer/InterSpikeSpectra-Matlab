function coeffs_new = debias_coefficients(coeffs, s, Theta)
% make a least-squares regression on the non-zero coefficients from the sparse regression
coeffs_new = coeffs;
non_zero_idx = find(coeffs);
coeffs_new(non_zero_idx) = Theta(:,non_zero_idx)\s;
end
