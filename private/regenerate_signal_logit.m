function reg = regenerate_signal_logit(Theta, coefs)
% regenerate a decomposed signal from basis functions and its coefficients

assert(size(Theta,2) == length(coefs),"Number of basis functions must match number of coefficients")
Xb  = exp(Theta*coefs);
reg = Xb./(1+Xb);

end