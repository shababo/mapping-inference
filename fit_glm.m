function fit_params = fit_glm(stims_x, stims_t, spikes)

num_basis_funcs = 10;
num_spatial_pos = 121;

num_params = 1 + 2*num_basis_funcs + num_spatial_pos;

bs = basisFactory.makeSmoothTemporalBasis('raised cosine', .750, num_basis_funcs, @(x) 1500);
basis = bs.B';

x0= rand(num_params,1);
obj_fun = @(x) glm_log_likelihood(x, stims_x, stims_t, spikes, basis);
options = optimoptions('fmincon','Algorithm','interior-point',...
                                'Display','iter','GradObj','on',...
                                'UseParallel',true,'Diagnostics','on');
                            
% options = optimoptions('fminunc','Algorithm','trust-region', 'Display','iter','GradObj','on');

% options = optimset('Display','iter');

ub = Inf*ones(size(x0));                            
lb = -Inf*ones(size(x0));
lb(2:num_spatial_pos+1) = 0;
                            

delete(gcp('nocreate'))
this_pool = parpool();
% fit_params = fminunc(obj_fun,x0,options)
fit_params = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
% fit_params = fminsearch(obj_fun,x0,options);
delete(this_pool)


function [log_likelihood, grad_ll] = glm_log_likelihood(param_vec, stims_x, stims_t, spikes, basis)

num_basis_funcs = size(basis,1);
num_trials = size(stims_t,1);
num_timebins = size(stims_t,2);

dt = 1/20000;

param_count = 1;
baseline_rate = param_vec(param_count); param_count = param_count + 1;

% build spatial layout of cell
num_spatial_pos = 121;
spatial_footprint = param_vec(param_count:param_count + num_spatial_pos - 1);
param_count = param_count + num_spatial_pos;

% a_stim = param_vec(param_count); param_count = param_count + 1;
% c_stim = param_vec(param_count); param_count = param_count + 1;

weights_stim = param_vec(param_count:param_count + num_basis_funcs - 1);
param_count = param_count + num_basis_funcs;
% for i = 1:num_basis_funcs
% %     phi = param_vec(param_count); param_count = param_count + 1;
%     weight =  param_vec(param_count); param_count = param_count + 1;
%     stim_filter = stim_filter + weight * basis(i,:); %raised_cosine(t,a_stim,c_stim,phi);
% end
stim_filter = weights_stim'*basis;

weights_hist = param_vec(param_count:param_count + num_basis_funcs - 1);
param_count = param_count + num_basis_funcs;
% a_hist = param_vec(param_count); param_count = param_count + 1;
% c_hist = param_vec(param_count); param_count = param_count + 1;
% for i = 1:num_basis_funcs
% %     phi = param_vec(param_count); param_count = param_count + 1;
%     weight =  param_vec(param_count); param_count = param_count + 1;
%     history_filter = history_filter + weight * basis(i,:); %raised_cosine(t,a_hist,c_hist,phi);
% end
history_filter = weights_hist'*basis;

total_filter_out = baseline_rate*ones(size(spikes));

% for grad
d_filter_out_baseline = baseline_rate*ones(size(spikes));
d_filter_out_spatial = zeros([size(spikes) num_spatial_pos]);
d_filter_out_stim_weights = zeros([size(spikes) num_basis_funcs]);
d_filter_out_hist_weights = zeros([size(spikes) num_basis_funcs]);

for i = 1:size(stims_t,1)
    
    stim_filter_conv_out = conv(stims_t(i,:),stim_filter);
    stim_filter_conv_out = stim_filter_conv_out(1:size(spikes,2));
    
    hist_filter_conv_out = conv(spikes(i,:),history_filter);
    hist_filter_conv_out = hist_filter_conv_out(1:size(spikes,2));
    
    total_filter_out(i,:) = total_filter_out(i,:) + ...
        spatial_footprint(stims_x(i,1)) * stims_x(i,2) * stim_filter_conv_out + ...
        hist_filter_conv_out;
        
    for j = 1:num_spatial_pos
        if stims_x(i,1) == j
            d_filter_out_spatial(i,:,j) = spatial_footprint(j) * ...
                stims_x(i,2) * stim_filter_conv_out;
        end
    end
    
    for j = 1:num_basis_funcs
        this_conv = conv(stims_t(i,:),weights_stim(j)*basis(j,:));
        this_conv = this_conv(1:num_timebins);
        d_filter_out_stim_weights(i,:,j) = ...
            spatial_footprint(stims_x(i,1)) * stims_x(i,2) * this_conv;
    end
    
    for j = 1:num_basis_funcs
        this_conv = conv(spikes(i,:),weights_hist(j)*basis(j,:));
        d_filter_out_hist_weights(i,:,j) = this_conv(1:num_timebins);
    end
    
end

lambda = dt*exp(total_filter_out);

log_likelihood = sum(sum(total_filter_out.*spikes - lambda));
log_likelihood = -log_likelihood; % because we are using fmin***

grad_ll = zeros(size(param_vec));

param_count = 1;

grad_ll(param_count) = sum(sum(-d_filter_out_baseline.*spikes + lambda.*d_filter_out_baseline));
param_count = param_count + 1;

for i = 1:num_spatial_pos
    grad_ll(param_count) = ...
        sum(sum(-d_filter_out_spatial(:,:,i).*spikes + lambda.*d_filter_out_spatial(:,:,i)));
    param_count = param_count + 1;
end

%stim filter
for i = 1:num_basis_funcs
    grad_ll(param_count) = ...
        sum(sum(-d_filter_out_stim_weights(:,:,i).*spikes + lambda.*d_filter_out_stim_weights(:,:,i)));
    param_count = param_count + 1;
end

%hist filter
for i = 1:num_basis_funcs
    grad_ll(param_count) = ...
        sum(sum(-d_filter_out_hist_weights(:,:,i).*spikes + lambda.*d_filter_out_hist_weights(:,:,i)));
    param_count = param_count + 1;
end




