function [spectrum, y] = compute_spectrum_according_to_actual_spectrum_lasso_logit(spectrum_i, s, Theta, lambda_max, alpha)
abs_tol = 1e-6;
max_iter = 10;  % precision
lambda_min = 0;
% bisection search
for i = 1:max_iter
    % try new lambda
    actual_lambda = lambda_min + (lambda_max - lambda_min)/2;  
    % make the regression with specific lambda
    y = lassoglm(Theta, s, 'binomial', 'Lambda', actual_lambda, 'Alpha', alpha);
    y = debias_coefficients(y, s, Theta); % coefficient debiasing
    spectrum = pool_frequencies(y, length(s));
    spectrum = spectrum ./ sum(spectrum);
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