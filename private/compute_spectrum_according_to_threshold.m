function [spectrum, rho] = compute_spectrum_according_to_threshold(s, Theta, threshold, tol, max_iter, verbose)
abs_tol = 1e-6;
% initial Lambda-step for estimating an upper bound for lambda
lambda_step = 0.5;
[lambda_max, lambda_min, y_act, rho_act] = find_lambda_max(s, Theta, lambda_step, threshold);
% bisection search
for i = 1:max_iter
    % check whether max iterations or tolerance-level reached
    if i == max_iter
        if rho_act > 1 - abs_tol
            spectrum_i = pool_frequencies(y_act, length(s));
            spectrum_i = spectrum_i ./ sum(spectrum_i);
            [spectrum, y] = compute_spectrum_according_to_actual_spectrum(spectrum_i, s, Theta, actual_lambda);
            rr = corrcoef(regenerate_signal(Theta, y), s);
            rho = rr(2);
            if verbose
                warning("Algorithm stopped due to maximum number of lambda's were tried without convergence to the specified `threshold`. Perfect Decomposition achieved.")
            end            
        else    
            spectrum = pool_frequencies(y_act, length(s));
            spectrum = spectrum ./ sum(spectrum);
            rho = rho_act;
            if verbose
                warning("Algorithm stopped due to maximum number of lambda's were tried without convergence to the specified `threshold`.")
            end
        end
        break
    elseif abs(rho_act - threshold) <= tol
        spectrum = pool_frequencies(y_act, length(s));
        spectrum = spectrum ./ sum(spectrum);
        rho = rho_act;
        break
    end
    % try new lambda
    actual_lambda = lambda_min + (lambda_max - lambda_min)/2;  
    % make the regression with specific lambda
    y = lasso(Theta, s, 'Lambda', actual_lambda);
    % check whether the regenerated signal matches with the given threshold
    rr = corrcoef(regenerate_signal(Theta, y), s);
    
    % pick the new bisection interval
    if isnan(rr(2))
        lambda_max = actual_lambda;
    elseif rr(2) < threshold
        lambda_max = actual_lambda;
        y_act(:) = y;
        rho_act = rr(2);
    elseif rr(2) > threshold
        lambda_min = actual_lambda;
        rho_act = rr(2);
        y_act(:) = y;
    end
end
end