function fit_params = fit_glm(stims_x, stims_t, spikes)

num_basis_funcs = 4;
num_spatial_pos = 121;

num_params = 1 + 2*(num_basis_funcs + 2) + num_spatial_pos;

problem.solver = 'fminunc';
problem.x0= zeros(num_params,1);
problem.objective = @(x) glm_log_likelihood(x, stims_x, stims_t, spikes);
problem.options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter');

fit_params = fminunc(problem);


function log_likelihood = glm_log_likelihood(param_vec, stims_x, stims_t, spikes)

num_basis_funcs = 4;

% build the filters from basis function
stim_filter = zeros(1,size(spikes,2));
history_filter = zeros(1,size(spikes,2));

t = (1:size(spikes,2))/20000;
param_count = 1;
baseline_rate = param_vec(param_count); param_count = param_count + 1;

% build spatial layout of cell
num_spatial_pos = 121;
sptial_footprint = param_vec(param_count:param_count + num_spatial_pos-1);
param_count = param_count + num_spatial_pos;

a_stim = param_vec(param_count); param_count = param_count + 1;
c_stim = param_vec(param_count); param_count = param_count + 1;

for i = 1:num_basis_funcs
    phi = param_vec(param_count); param_count = param_count + 1;
    stim_filter = stim_filter + raised_cosine(t,a_stim,c_stim,phi);
end

a_hist = param_vec(param_count); param_count = param_count + 1;
c_hist = param_vec(param_count); param_count = param_count + 1;
for i = 1:num_basis_funcs
    phi = param_vec(param_count); param_count = param_count + 1;
    history_filter = history_filter + raised_cosine(t,a_hist,c_hist,phi);
end

linear_filter_out = baseline_rate*ones(size(spikes));

for i = 1:size(stims_t,1)
    
    stim_filter_conv_out = conv(stims_t(i,:),stim_filter);
    stim_filter_conv_out = stim_filter_conv_out(1:size(stims_t,2));
    stim_filter_conv_out = stim_filter_conv_out(1:size(spikes,2));
    
    hist_filter_conv_out = conv(spikes(i,:),history_filter);
    hist_filter_conv_out = hist_filter_conv_out(1:size(stims_t,2));
    hist_filter_conv_out = hist_filter_conv_out(1:size(spikes,2));
    
    linear_filter_out(i,:) = linear_filter_out(i,:) + ...
        sptial_footprint(stims_x(i,1)) * stims_x(i,2) * stim_filter_conv_out + ...
        hist_filter_conv_out;
    
end

lambda = exp(linear_filter_out);

log_likelihood = sum(sum(linear_filter_out.*spikes)) - (1/20000)*sum(sum(lambda));
log_likelihood = -log_likelihood; % because we are using fmin***





