function out = raised_cosine(x,a,c,phi)

shift_term = a * log(x + c);
indicator_vec = phi - pi <= shift_term & shift_term <= phi + pi;

out = 0.5 * cos(a * log(x + c) - phi) + 0.5;
out = out .* indicator_vec;



