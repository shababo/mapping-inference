function spikes = sim_glm(params, stims_x, stims_t, dt)

spikes = zeros(size(stims_t));
total_filter_out = params.baseline * ones(size(stims_t));

for i = 1:size(stims_t,1)
    
    stim_filter_conv_out = conv(stims_t(i,:),params.stim_filter);
    stim_filter_conv_out = stim_filter_conv_out(1:size(stims_t,2));
    
    total_filter_out(i,:) = total_filter_out(i,:) + ...
        params.spatial_footprint(stims_x(i,1)) * stim_filter_conv_out + ...
        hist_filter_conv_out;
    
end

for i = 1:sizse(stims_t,2)
    

lambda = dt*exp(total_filter_out);

spikes = poissrnd(lambda);

