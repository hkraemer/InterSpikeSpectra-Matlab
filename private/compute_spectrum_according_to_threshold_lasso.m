function [spectrum, rho] = compute_spectrum_according_to_threshold_lasso(s, Theta, threshold, tol, max_iter, alpha, verbose)
abs_tol = 1e-6;
% initial Lambda-step for estimating an upper bound for lambda
lambda_step = 0.5;
[lambda_max, lambda_min, y_act, ~] = find_lambda_max(s, Theta, lambda_step, threshold, alpha);

% bisection search
for i = 1:max_iter    
    % try new lambda
    actual_lambda = lambda_min + (lambda_max - lambda_min)/2;  
    % make the regression with specific lambda
    y = lasso(Theta, s, 'Lambda', actual_lambda, 'Alpha', alpha);
    % check whether the regenerated signal matches with the given threshold
    rr = corrcoef(regenerate_signal(Theta, y), s);
    
    % pick the new bisection interval
    if isnan(rr(2))
        lambda_max = actual_lambda;
    elseif rr(2) < threshold
        lambda_max = actual_lambda;
        y_act(:) = y;
    elseif rr(2) > threshold
        lambda_min = actual_lambda;
        y_act(:) = y;
    end
    rho_act = rr(2);
    % check whether max iterations or tolerance-level reached
    if i == max_iter
        if rho_act > 1 - abs_tol
            if verbose
                warning("Algorithm stopped due to maximum number of lambda's were tried without convergence to the specified `threshold`. Perfect Decomposition achieved.")
            end 
            y_act = debias_coefficients(y_act, s, Theta); % coefficient debiasing
            spectrum_i = pool_frequencies(y_act, length(s));
            spectrum_i = spectrum_i ./ sum(spectrum_i); % normalization
            [spectrum, y] = compute_spectrum_according_to_actual_spectrum_lasso(spectrum_i, s, Theta, actual_lambda, alpha);
            rr = corrcoef(regenerate_signal(Theta, y), s);
            rho = rr(2); 
        else
            if verbose
                warning("Algorithm stopped due to maximum number of lambda's were tried without convergence to the specified `threshold`.")
            end
            if isnan(rho_act) % account for case where the last iterations yielded a NaN-value
                y_act(:) = lasso(Theta, s, 'Lambda', lambda_min, 'Alpha', alpha);
                % check whether the regenerated signal matches with the given threshold
                rr = corrcoef(regenerate_signal_logit(Theta, y_act), s);
                rho_act = rr(2);
            end 
            y_act = debias_coefficients(y_act, s, Theta); % coefficient debiasing
            spectrum = pool_frequencies(y_act, length(s));
            spectrum = spectrum ./ sum(spectrum); % normalization
            rho = rho_act;
        end
        break
    elseif abs(rho_act - threshold) <= tol
        y_act = debias_coefficients(y_act, s, Theta); % coefficient debiasing
        spectrum = pool_frequencies(y_act, length(s));
        spectrum = spectrum ./ sum(spectrum); % normalization
        rho = rho_act;        
        break
    end
end
end