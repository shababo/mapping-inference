function [prior_params]=initialize_prior(prior_path)
% specify one set of parameters

sum_path = [prior_path];
load(sum_path)

bounds=struct;
bounds.PR=[0.01 1];
bounds.gain=[0.005 0.06];
bounds.delay_mu=[0 60];
bounds.delay_sigma=[0.1 10];





ini_mean=0;
ini_sigma=2;
prior_params=struct;

prior_params.GP_params=summary_results;
prior_params.boundary_params= [30 30 70];

prior_params.shift_x.dist='normal';
prior_params.shift_x.type='individual';
prior_params.shift_x.mean=summary_results.x.shift_params.mean;
prior_params.shift_x.log_sigma=log(summary_results.x.shift_params.var)/2;

prior_params.shift_y.dist='normal';
prior_params.shift_y.type='individual';
prior_params.shift_y.mean=summary_results.y.shift_params.mean;
prior_params.shift_y.log_sigma=log(summary_results.y.shift_params.var)/2;

prior_params.shift_z.dist='normal';
prior_params.shift_z.type='individual';
prior_params.shift_z.mean=summary_results.z.shift_params.mean;
prior_params.shift_z.log_sigma=log(summary_results.z.shift_params.var)/2;

prior_params.gain.dist='logit-normal';
prior_params.gain.type='individual';
prior_params.gain.mean=ini_mean;
prior_params.gain.log_sigma=log(ini_sigma);
prior_params.gain.bounds.up=bounds.gain(2);
prior_params.gain.bounds.low=bounds.gain(1);



prior_params.PR.dist='logit-normal';
prior_params.PR.type='individual';
prior_params.PR.mean=ini_mean;
prior_params.PR.log_sigma=log(ini_sigma);
prior_params.PR.bounds.up=bounds.PR(2);
prior_params.PR.bounds.low=bounds.PR(1);
%         variational_params.PR.bounds.prob_logit=;

prior_params.delay_mu.dist='logit-normal';
prior_params.delay_mu.type='individual';
prior_params.delay_mu.mean=ini_mean;
prior_params.delay_mu.log_sigma=log(ini_sigma);
prior_params.delay_mu.bounds.up=bounds.delay_mu(2);
prior_params.delay_mu.bounds.low=bounds.delay_mu(1);


prior_params.delay_sigma.dist='logit-normal';
prior_params.delay_sigma.type='individual';
prior_params.delay_sigma.mean=ini_mean;
prior_params.delay_sigma.log_sigma=log(ini_sigma);
prior_params.delay_sigma.bounds.up=bounds.delay_sigma(2);
prior_params.delay_sigma.bounds.low=bounds.delay_sigma(1);


prior_params.shapes=struct;
prior_params.shapes.dist='normal';
prior_params.shapes.type='individual';
prior_params.shapes.locations=zeros(0,3);
prior_params.shapes.mean=zeros(0,1);
prior_params.shapes.log_sigma=zeros(0,1);
