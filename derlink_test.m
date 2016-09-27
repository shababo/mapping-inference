function dmu = derlink_test(mu)

linear_i = mu >= 1;
exp_i = mu < 1;

dmu = zeros(size(mu));
dmu(linear_i) = 1.0;
dmu(exp_i) = 1./mu(exp_i);