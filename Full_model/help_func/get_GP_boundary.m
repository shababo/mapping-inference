function [pred_mean, post_mean] = get_GP_boundary(Xstar, params)
X=params.X;
Y=params.Y;
boundary=params.boundary;
buffer=params.buffer;
mag=params.mag;
tau=params.tau;
epsilon=1e-5;
noninfo_var_func = ...
    @(x) (mag.*((abs(x)<(boundary-buffer) ))  + mag.*exp( -sqrt(abs(x)-boundary+buffer)).*(abs(x)>(boundary-buffer-epsilon) ));
%% Get the posterior mean:
sigma_var= noninfo_var_func(X);

nsq=sum(X.^2,2);
K=bsxfun(@plus,nsq,nsq');
K=bsxfun(@minus,K,(2*X)*X.');
sigma_mat = sigma_var*ones(1,length(X));
K=sigma_mat.*exp(-K/tau).*sigma_mat';
K_exp=K+ diag(ones(length(X),1))*0.05;

% Correlation between new points and the old points:
nsqstar=sum(Xstar.^2,2);
K_pred = nsqstar*ones(1,length(nsq)) +  ones(length(nsqstar),1)*nsq';
K_pred=bsxfun(@minus,K_pred,(2*Xstar)*X.');

sigma_varstar= noninfo_var_func(Xstar);
sigma_mat_star_left= sigma_varstar*ones(1,length(nsq));
sigma_mat_star_right= ones(length(nsqstar),1)*sigma_var';

K_pred= sigma_mat_star_left.*exp(-K_pred/tau).*sigma_mat_star_right;


post_mean= K*inv(K_exp)*(Y);
pred_mean= K_pred*inv(K_exp)*(Y);

