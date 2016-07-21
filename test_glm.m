% load 'work_spiking_data_0720.mat'


stims_x = [repmat((1:121)',5*3,1) [25*ones(121*5,1); 50*ones(121*5,1); 100*ones(121*5,1);]];

stim_t = zeros(1,1500);
stim_t(100) = 1;
stims_t = repmat(stim_t,size(stims_x,1),1);

spikes = zeros(size(stims_t));
count = 1;
for i = 1:4:6
    for j = 1:11
        for k = 1:11
            for l = 1:5
                spike_inds = find(detection_grids_3_31_s2c2_r2_3{i}{j,k}(l,1:1500));
                spikes(count,spike_inds) = 1;
                count = count + 1;
            end
        end
    end
end

%%

fit_params = fit_glm(stims_x, stims_t, spikes);

%%
num_basis_funcs = 4;

% build the filters from basis function
stim_filter = zeros(1,size(spikes,2));
history_filter = zeros(1,size(spikes,2));

param_count = 1;
baseline_rate = fit_params(param_count); param_count = param_count + 1;

t = (1:1500)/20000;

% build spatial layout of cell
num_spatial_pos = 121;
sptial_footprint = fit_params(param_count:param_count + num_spatial_pos-1);
param_count = param_count + num_spatial_pos;

a_stim = fit_params(param_count); param_count = param_count + 1;
c_stim = fit_params(param_count); param_count = param_count + 1;

for i = 1:num_basis_funcs
    phi = fit_params(param_count); param_count = param_count + 1;
    stim_filter = stim_filter + raised_cosine(t,a_stim,c_stim,phi);
end

a_hist = fit_params(param_count); param_count = param_count + 1;
c_hist = fit_params(param_count); param_count = param_count + 1;
for i = 1:num_basis_funcs
    phi = fit_params(param_count); param_count = param_count + 1;
    history_filter = history_filter + raised_cosine(t,a_hist,c_hist,phi);
end

figure
subplot(311)
imagesc(reshape(sptial_footprint,11,11))

subplot(312)
plot(stim_filter)

subplot(313)
plot(history_filter)
