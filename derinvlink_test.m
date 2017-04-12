function out = derinvlink_test(resp)

linear_i = resp >= 1;
exp_i = resp < 1;

out = zeros(size(resp));
out(linear_i) = 1.0;
out(exp_i) = resp(exp_i);