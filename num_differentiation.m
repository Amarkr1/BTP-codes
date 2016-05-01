function [value] = num_differentiation(w,h,lambda,String)
delta_lambda = 1e-3;
val = calc_neff(w,h,lambda+delta_lambda,String) - calc_neff(w,h,lambda-delta_lambda,String);
value = val/(2*delta_lambda);
end
