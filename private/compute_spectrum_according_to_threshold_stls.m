function [spectrum, rho] = compute_spectrum_according_to_threshold_stls(s, Theta, threshold, tol, max_iter, verbose)
% Determine the right STLS-regularization parameter with respect to ρ_thres & tol using Bisection-search.
abs_tol = 1e-6;
% lambda_max/ lambda_min correspond to the maximum/minimum value in the coeffs of least squares
y_act = pinv(Theta)*s;    % least squares
%y_act = Theta\s;
lambda_maxx = max(y_act);
lambda_max = max(y_act);
lambda_min = min(y_act);

% bisection search
for i = 1:max_iter    
    % try new lambda
    actual_lambda = lambda_min + (lambda_max - lambda_min)/2;  
    % make the regression with specific lambda
    y_act(:) = stls(s, Theta, actual_lambda);
    % check whether the regenerated signal matches with the given threshold
    rr = corrcoef(regenerate_signal(Theta, y_act), s);
    
    % pick the new bisection interval
    if isnan(rr(2))
        lambda_max = actual_lambda;
    elseif rr(2) < threshold
        lambda_max = actual_lambda;
    elseif rr(2) > threshold
        lambda_min = actual_lambda;
    end
    rho_act = rr(2);
    % check whether max iterations or tolerance-level reached
    if i == max_iter
        if rho_act > 1 - abs_tol
            if verbose
                warning("Algorithm stopped due to maximum number of lambda's were tried without convergence to the specified `threshold`. Perfect Decomposition achieved.")
            end   
            spectrum_i = pool_frequencies(y_act, length(s));
            spectrum_i = spectrum_i ./ sum(spectrum_i); % normalization
            [spectrum, y] = compute_spectrum_according_to_actual_spectrum_stls(spectrum_i, s, Theta, actual_lambda, lambda_maxx);
            rr = corrcoef(regenerate_signal(Theta, y), s);
            rho = rr(2);
        else
            if verbose
                warning("Algorithm stopped due to maximum number of lambda's were tried without convergence to the specified `threshold`.")
            end
            if isnan(rho_act) % account for case where the last iterations yielded a NaN-value
                y_act(:) = stls(s, Theta, lambda_min);
                % check whether the regenerated signal matches with the given threshold
                rr = corrcoef(regenerate_signal(Theta, y_act), s);
                rho_act = rr(2);
            end
            spectrum = pool_frequencies(y_act, length(s));
            spectrum = spectrum ./ sum(spectrum); % normalization
            rho = rho_act;
        end
        break
    elseif abs(rho_act - threshold) <= tol
        spectrum = pool_frequencies(y_act, length(s));
        spectrum = spectrum ./ sum(spectrum); % normalization
        rho = rho_act;
        break
    end
end
end