function [pred_mean, post_mean,prior_var] = get_2DGP_boundary(Xstar, params)
%%
X=params.X;
Y=params.Y;
epsilon=1e-5;
sigma_noise=params.x.sigma_noise;
x_boundary=params.x.boundary;
x_buffer=params.x.buffer;
x_mag=params.x.mag;
x_tau=params.x.tau;
x_noninfo_var_func = @(x) x_mag*( (1./(1+exp(-x_boundary-x ))).*(1./(1+exp(x-x_boundary))));
y_boundary=params.y.boundary;
y_buffer=params.y.buffer;
y_mag=params.y.mag;
y_tau=params.y.tau;
y_noninfo_var_func = @(x) y_mag*( (1./(1+exp(-y_boundary-x ))).*(1./(1+exp(x-y_boundary))));

distance_factor=zeros(1,2);
distance_factor(1)=params.x.distance_factor;
distance_factor(2)=params.x.distance_factor;

tau_factor=zeros(1,2);
tau_factor(1)=params.x.tau_factor;
tau_factor(2)=params.y.tau_factor;

%% Shrink the X and Xstar
boundary=zeros(1,2);
boundary(1)=    params.x.shrink_bound;
boundary(2)=    params.y.shrink_bound;

for i = 1:2
    tail_indices=abs(X(:,i))>boundary(i);
    tmp=X(tail_indices,i);
    mag_tmp=(tau_factor(i))*(abs(tmp)-boundary(i)).^distance_factor(i) +boundary(i);
    X(tail_indices,i)=sign(tmp).*mag_tmp;
    
    tail_indices=abs(Xstar(:,i))>boundary(i);
    tmp=Xstar(tail_indices,i);
    mag_tmp=(tau_factor(i))*(abs(tmp)-boundary(i)).^distance_factor(i) +boundary(i);
    Xstar(tail_indices,i)=sign(tmp).*mag_tmp;
end

%% Get the posterior mean:
sigma_var(:,1)= max(epsilon,x_noninfo_var_func(X(:,1)));
sigma_var(:,2)= max(epsilon,y_noninfo_var_func(X(:,2)));
taus=[x_tau y_tau];
K_exp=cell([2 1]);K_data=cell([2 1]);
for i=1:2
    this_X=X(:,i);
    sigma_mat = sigma_var(:,i)*ones(1,length(this_X));
    K_data{i}=sigma_mat.*get_kernel_cor(this_X,this_X,taus(i)).*sigma_mat';   
end

% Correlation between new points and the old points:
K_star=cell([2 1]);
clear('sigma_varstar')
sigma_varstar(:,1)= max(epsilon,x_noninfo_var_func(Xstar(:,1)));
sigma_varstar(:,2)= max(epsilon,y_noninfo_var_func(Xstar(:,2)));    
for i=1:2
    this_X=X(:,i);
    this_Xstar=Xstar(:,i);
    sigma_mat_star_left= sigma_varstar(:,i)*ones(1,length(this_X));
    sigma_mat_star_right= ones(length(this_Xstar),1)*sigma_var(:,i)';
    K_star{i}= sigma_mat_star_left.*get_kernel_cor(this_Xstar,this_X,taus(i)).*sigma_mat_star_right;
end

K_star_xy=K_star{1}.*K_star{2};
K_data_xy=K_data{1}.*K_data{2};
K_exp_xy=K_data_xy+ diag(ones(length(this_X),1))*sigma_noise;

K_inv=inv(K_exp_xy);


prior_mean=params.prior_mean*ones(length(Y),1);
KY=K_inv*(Y-prior_mean);


post_mean=prior_mean+K_data_xy*KY;
pred_mean=params.prior_mean*ones(length(Xstar),1)+K_star_xy*KY;


prior_var = sigma_varstar(:,1).*sigma_varstar(:,2);
