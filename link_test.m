function resp = link_test(mu)

linear_i = mu >= 1;
exp_i = mu < 1;

resp = zeros(size(mu));
resp(linear_i) = (-1 + mu(linear_i));
resp(exp_i) = log(mu(exp_i));