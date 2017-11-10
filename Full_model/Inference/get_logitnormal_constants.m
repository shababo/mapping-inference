function [logitnormal_constants] = get_logitnormal_constants(alpha,beta)
beta=exp(beta);
logitnormal_constants=struct;
logitnormal_constants.mu_by_sigma2 =[alpha./(beta.^2)]';
logitnormal_constants.sigma_inv= [1./(beta)]';
logitnormal_constants.sigma_inv2= [1./(beta.^2)]';
logitnormal_constants.sigma_inv3= [1./(beta.^3)]';
end




