
function [log_likelihood, grad_ll] = glm_log_likelihood(param_vec, stims_x, stims_t, spikes, basis)

% global g_stims_x
% global g_stims_t
% global g_spikes
% global g_basis
% 
% stims_x = g_stims_x;
% stims_t = g_stims_t;
% spikes = g_spikes;
% basis = g_basis;

num_basis_funcs = size(basis,1);
num_trials = size(stims_t,1);
num_timebins = size(stims_t,2);

dt = 1/1000;

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

total_filter_out = baseline_rate + zeros(size(spikes));

% for grad
d_filter_out_baseline = ones(size(spikes));
d_filter_out_spatial = zeros([size(spikes) num_spatial_pos]);
d_filter_out_stim_weights = zeros([size(spikes) num_basis_funcs]);
d_filter_out_hist_weights = zeros([size(spikes) num_basis_funcs]);

for i = 1:size(stims_t,1)
    
%     t = convmtx(stims_t(i,:),length(stim_filter),length(stim_filter));
    stim_filter_conv_out = conv(stims_t(i,:),stim_filter);
    stim_filter_conv_out = stim_filter_conv_out(1:size(spikes,2));
    
%     if i == 1
%         figure(123)
%         plot(stim_filter_conv_out)
%         figure(124)
%         plot(stim_filter)
%         drawnow
%     end
    
    hist_filter_conv_out = conv(spikes(i,:),history_filter);
    hist_filter_conv_out = hist_filter_conv_out(1:size(spikes,2));
    
%     
%     figure(124)
%     plot(hist_filter_conv_out)
    
    total_filter_out(i,:) = total_filter_out(i,:) + ...
        spatial_footprint(stims_x(i,1)) * stim_filter_conv_out + ...
        hist_filter_conv_out;
        
    for j = 1:num_spatial_pos
        if stims_x(i,1) == j
            d_filter_out_spatial(i,:,j) = ...
                stim_filter_conv_out;
        end
    end
    
    for j = 1:num_basis_funcs
        this_conv = conv(stims_t(i,:),basis(j,:));
        this_conv = this_conv(1:num_timebins);
        d_filter_out_stim_weights(i,:,j) = ...
            spatial_footprint(stims_x(i,1)) * this_conv;
    end
    
    for j = 1:num_basis_funcs
        this_conv = conv(spikes(i,:),basis(j,:));
        d_filter_out_hist_weights(i,:,j) = this_conv(1:num_timebins);
    end
    
end

lambda = dt*exp(total_filter_out);



log_likelihood = sum(sum(-total_filter_out.*spikes + lambda));

% if isinf(log_likelihood)
%     log_likelihood = 1e200;
% end

grad_ll = zeros(size(param_vec));

% d_filter_out_baseline = d_filter_out_baseline/numel(d_filter_out_stim_weights);
% d_filter_out_spatial = d_filter_out_spatial/numel(d_filter_out_stim_weights);
% d_filter_out_stim_weights = d_filter_out_stim_weights/numel(d_filter_out_stim_weights);
% d_filter_out_hist_weights = d_filter_out_hist_weights/numel(d_filter_out_stim_weights);

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

% grad_ll = grad_ll ;

% log_likelihood
% 
% if any(isnan(grad_ll))
%     disp('g NaN')
% end
% if any(isinf(grad_ll))
%     disp('g Inf')
% end

% assignin('base','total_filter_out',total_filter_out);
% assignin('base','lambda',lambda);
% assignin('base','stim_filter_conv_out',stim_filter_conv_out);
% assignin('base','hist_filter_conv_out',hist_filter_conv_out);
% assignin('base','d_filter_out_spatial',d_filter_out_spatial);
% assignin('base','d_filter_out_stim_weights',d_filter_out_stim_weights);
% assignin('base','d_filter_out_hist_weights',d_filter_out_hist_weights);
