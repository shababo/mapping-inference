function [figure_handle]=shape_same_plot(shifts, neurons, figure_handle,method,tau, sigma_gp, sigma,...
    mean_func)


n_cell=length(neurons);

corrected_grid=cell([n_cell 1]);
for i_cell = 1:n_cell
    corrected_grid{i_cell}=[neurons(i_cell).stim_grid - shifts(i_cell)];
end

% plot(corrected_grid], );

X=[corrected_grid{:}]';
Y=[neurons(:).scaled_current]';
nsq=sum(X.^2,2);
K=bsxfun(@plus,nsq,nsq');
K=bsxfun(@minus,K,(2*X)*X.');
K=sigma_gp*exp(-K/tau);

K_exp=K+ diag(ones(length([corrected_grid{:}]),1))*sigma;

post_mean= K*inv(K_exp)*(Y-mean_func(X))+mean_func(X);
post_cov=K-K*inv(K_exp)*K;

% scatter(X,post_mean)
[uniq_grid, ix]=sort(X);
post_var=diag(post_cov);

plot(uniq_grid,post_mean(ix))
hold on;
for i_cell=1:n_cell
    scatter(neurons(i_cell).stim_grid- shifts(i_cell),neurons(i_cell).scaled_current);
    hold on;
end
%plot(uniq_grid,sqrt(post_cov(ix)+sigma^2))
% Draw confidence band?

xlim([-70 70])

xlabel('Stimulation location (adjusted)')
ylabel('Scaled current')
title([method ' ' 'Tau ' num2str(round(tau,1)) ', ' 'Sigma ' num2str(round(sigma_gp,3))])

