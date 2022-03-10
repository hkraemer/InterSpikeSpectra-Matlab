function [spectrum, y_act] = compute_spectrum_according_to_actual_spectrum_stls(spectrum_i, s, Theta, lambda_min, lambda_max)
assert(lambda_min<lambda_max)

abs_tol = 1e-6;
max_iter = 10;  % precision
y_act = zeros(size(Theta,2),1);     % preallocation
% bisection search
for i = 1:max_iter
    % try new lambda
    actual_lambda = lambda_min + (lambda_max - lambda_min)/2;  
    % make the regression with specific lambda
    y_act(:) = stls(s, Theta, actual_lambda);
    spectrum = pool_frequencies(y_act, length(s));
    spectrum = spectrum ./ sum(spectrum); % normalization
    % check whether the spectrum matches with the initial spectrum (input)
    rr = corrcoef(spectrum, spectrum_i);
    
    % pick the new bisection interval
    if rr(2) > 1-abs_tol
        lambda_max = actual_lambda;
    elseif rr(2) < 1-abs_tol
        lambda_min = actual_lambda;
    end
end
end