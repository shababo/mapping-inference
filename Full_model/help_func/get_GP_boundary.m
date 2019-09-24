function [pred_mean, post_mean,prior_var] = get_GP_boundary(Xstar, params)
%%
X=params.X;
Y=params.Y;
boundary=params.boundary;
buffer=params.buffer;
mag=params.mag;
tau=params.tau;
distance_factor=params.distance_factor;

tau_factor=params.tau_factor;
epsilon=1e-5;
sigma_noise=params.sigma_noise;
noninfo_var_func = @(x) mag*( (1./(1+exp(-boundary-x ))).*(1./(1+exp(x-boundary))));
%% Shrink the X and Xstar
boundary=    params.shrink_bound;
tail_indices=abs(X)>boundary;
tmp=X(tail_indices);
 mag_tmp=(tau_factor)*(abs(tmp)-boundary).^distance_factor +boundary;
 X(tail_indices)=sign(tmp).*mag_tmp;
 
tail_indices=abs(Xstar)>boundary;
tmp=Xstar(tail_indices);
 mag_tmp=(tau_factor)*(abs(tmp)-boundary).^distance_factor +boundary;
 Xstar(tail_indices)=sign(tmp).*mag_tmp;
 
%% Get the posterior mean:
sigma_var= max(epsilon,noninfo_var_func(X));
% X=X./(epsilon+abs(X).^(distance_factor));
nsq=sum(X.^2,2);
K=bsxfun(@plus,nsq,nsq');
K=bsxfun(@minus,K,(2*X)*X.');
sigma_mat = sigma_var*ones(1,length(X));
% tau_mat = ( (bsxfun(@plus,nsq,nsq')).^(distance_factor)+1)*tau;
% A1=nsq*ones(1,length(nsq)); A2= ones(length(nsq),1)*nsq';
% tau_mat = (((abs(A1+A2)-abs(A1-A2))/2).^(distance_factor)+1)*tau;
K=sigma_mat.*exp(-K./tau).*sigma_mat';
K_exp=K+ diag(ones(length(X),1))*sigma_noise;

% Correlation between new points and the old points:
% Xstar=Xstar./(abs(Xstar).^(distance_factor)+epsilon);
nsqstar=sum(Xstar.^2,2);
K_pred = nsqstar*ones(1,length(nsq)) +  ones(length(nsqstar),1)*nsq';
K_pred=bsxfun(@minus,K_pred,(2*Xstar)*X.');
sigma_varstar= max(epsilon,noninfo_var_func(Xstar));
sigma_mat_star_left= sigma_varstar*ones(1,length(nsq));
sigma_mat_star_right= ones(length(nsqstar),1)*sigma_var';
% Find the minimun:
% A1=nsqstar*ones(1,length(nsq)) ; A2= ones(length(nsqstar),1)*nsq';
% tau_mat =( ((abs(A1+A2)-abs(A1-A2))/2).^(distance_factor)+1)*tau;
K_pred= sigma_mat_star_left.*exp(-K_pred./tau).*sigma_mat_star_right;

% With non-zero means: 
prior_mean=params.prior_mean*ones(length(Y),1);
post_mean=prior_mean+K*inv(K_exp)*(Y-prior_mean);
pred_mean= params.prior_mean*ones(length(Xstar),1)+K_pred*inv(K_exp)*(Y-prior_mean);
prior_var = sigma_varstar;
