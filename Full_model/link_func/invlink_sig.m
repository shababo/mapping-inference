function mu = invlink_sig(resp)
mu =(1-0.01)*sigmf(resp,[3 -1]);

% %%
% figure(1)
% plot(exp(nvec))
% hold on;
%  plot(sigmf(nvec,[3 -1]) )