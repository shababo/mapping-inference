function mu = invlink_test(resp)

linear_i = resp >= 0;
exp_i = resp < 0;

mu = zeros(size(resp));
mu(linear_i) = (1 + resp(linear_i));
mu(exp_i) = exp(resp(exp_i));