% load 'work_spiking_data_0720.mat'



stims_x = [repmat((1:121)',5*3,1) [25*ones(121*5,1); 50*ones(121*5,1); 100*ones(121*5,1);]];
num_spatial_pos = 121;
num_trials = 121*5*3;

num_cells = 10;

stims_x_value = repmat((1:121)',5*3,1);
stims_x_vec = zeros(num_trials,121);

powers = [25 50 100];

stims_t_downres = zeros(num_trials,75);

count = 1;
for i = 1:3
    for l = 1:5
        for j = 1:11
            for k = 1:11
                stims_t_downres(count,5:9) = powers(i);
                stims_x_vec(count,stims_x_value(count)) = 1;
                count = count + 1;
            end
        end
    end
end
%%

spikes_downres = cell(num_cells,1);
spatial_inits = zeros(11,11,num_cells);


for m = 1:num_cells
    

    % stim_t = zeros(size(stims_x,1),1500);
    % % stim_t(100:199) = 1;
    % % stims_t = repmat(stim_t,size(stims_x,1),1);

    spikes_downres{m} = zeros(size(stims_t_downres));

    count = 1;
    for i = 1:3
        for l = 1:5
            for j = 1:11
                for k = 1:11
                    if l <= size(all_detection_grids{m}{i}{j,k},1)
                        spike_inds = ceil(find(all_detection_grids{m}{i}{j,k}(l,1:1500))/20);
                        spikes_downres{m}(count,spike_inds) = 1;
                        if ~isempty(spike_inds)
                            spatial_inits(j,k,m) = spatial_inits(j,k,m) + 1;
                        end
                    end
                    count = count + 1;
                end 
            end
        end
    end

end

%%
% num_cells = 10;
% fit_params_tmp = fit_params;
clc

delete(gcp('nocreate'))
this_pool = parpool();
% % init_vals = fit_params;
fit_params_all_cells3 = cell(num_cells,1);
% fit_params{1} = init_vals;
for m = [1 3:10]
%     try
        init_values = 1e-9*ones(142,1);
        spatial_init_cell = spatial_inits(:,:,m)';
        spatial_init_cell = spatial_init_cell(:)/max(spatial_init_cell(:));
        init_values(2:num_spatial_pos+1) = init_values(2:num_spatial_pos+1) + spatial_init_cell*5e-3;
        if m > 1
            init_values(123:132) = fit_params_all_cells3{1}(123:132);
        end
        fit_params_all_cells3{m} = fit_glm(stims_x, stims_t_downres, spikes_downres{m}, init_values);
%     catch e
%         disp(['iter' num2str(m) 'failed'])
%     end
end

delete(this_pool)
%%

for i = [1]
    
    fit_params = fit_params_all_cells3{i};
    
    num_basis_funcs = 10;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', .750, num_basis_funcs, @(x) 75);
    basis = bs.B(1:end,:)';

    % build the filters from basis function
    stim_filter = zeros(1,size(spikes,2));
    history_filter = zeros(1,size(spikes,2));

    param_count = 1;
    baseline_rate = fit_params(param_count); param_count = param_count + 1;

    t = (0:74)/1000;

    % build spatial layout of cell
    num_spatial_pos = 121;
    sptial_footprint = fit_params(param_count:param_count + num_spatial_pos-1);
    param_count = param_count + num_spatial_pos;

    % a_stim = fit_params(param_count); param_count = param_count + 1;
    % c_stim = fit_params(param_count); param_count = param_count + 1;

    for j = 1:num_basis_funcs
    %     phi = fit_params(param_count); param_count = param_count + 1;
        weight =  fit_params(param_count); param_count = param_count + 1;
        stim_filter = stim_filter + weight * basis(j,:); %raised_cosine(t,a_stim,c_stim,phi);
    end

    % a_hist = fit_params(param_count); param_count = param_count + 1;
    % c_hist = fit_params(param_count); param_count = param_count + 1;
    for j = 1:num_basis_funcs
    %     phi = fit_params(param_count); param_count = param_count + 1;
        weight =  fit_params(param_count); param_count = param_count + 1;
        history_filter = history_filter + weight * basis(j,:); %raised_cosine(t,a_hist,c_hist,phi);
    end

    figure
    subplot(2,2,[1 3])
    imagesc(reshape(sptial_footprint,11,11)')
    title(['Cell ' num2str(i) ': Spatial Stim Filter'])
    axis off

    subplot(222)
    plot(t,stim_filter)
    title('Temporal Stim Filter')
    set(gca,'YTickLabels',{})

    subplot(224)
    plot(t,history_filter)
    title('Spike History/Intrinsic Firing Property Filter')
    set(gca,'YTickLabels',{})
    
    set(gcf,'Position',[1 1 1231 433])
    set(gcf, 'Color', 'w');
    export_fig(gcf, ['cell' num2str(i) '_glm_results_01.png'],'-png')
    

end

%% sim glm

for i = 1
    
    fit_params = fit_params_all_cells3{i};
    
    num_basis_funcs = 10;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', .750, num_basis_funcs, @(x) 75);
    basis = bs.B(1:end,:)';

    % build the filters from basis function
    params.stim_filter = zeros(1,size(spikes,2));
    params.hist_filter = zeros(1,size(spikes,2));

    param_count = 1;
    params.baseline = fit_params(param_count); param_count = param_count + 1;

    t = (0:74)/1000;

    % build spatial layout of cell
    num_spatial_pos = 121;
    params.spatial_footprint = fit_params(param_count:param_count + num_spatial_pos-1);
    param_count = param_count + num_spatial_pos;

    % a_stim = fit_params(param_count); param_count = param_count + 1;
    % c_stim = fit_params(param_count); param_count = param_count + 1;

    for j = 1:num_basis_funcs
    %     phi = fit_params(param_count); param_count = param_count + 1;
        weight =  fit_params(param_count); param_count = param_count + 1;
        params.stim_filter = params.stim_filter + weight * basis(j,:); %raised_cosine(t,a_stim,c_stim,phi);
    end
    figure; plot(params.stim_filter)
    
    % a_hist = fit_params(param_count); param_count = param_count + 1;
    % c_hist = fit_params(param_count); param_count = param_count + 1;
    for j = 1:num_basis_funcs
    %     phi = fit_params(param_count); param_count = param_count + 1;
        weight =  fit_params(param_count); param_count = param_count + 1;
        params.hist_filter = params.hist_filter + weight * basis(j,:); %raised_cosine(t,a_hist,c_hist,phi);
    end
    figure; plot(params.hist_filter)
    
    spikes_sim = sim_glm(params,stims_x,stims_t_downres,1/1000);
    
    
    spikes_grids = cell(3,1);
    
    count = 1;
    for ii = 1:length(spike_grids)
        spikes_grids{ii} = cell(11,11);
        for ll = 1:5
            for jj = 1:11
                for kk = 1:11
                    if jj == 1 && kk == 1
                        spikes_grids{ii}{jj,kk} = zeros(5,size(stims_t_downres,2));
                    end
                    spikes_grids{ii}{jj,kk}(ll,:) = spikes_sim(count,:)*70;
                    count = count + 1;
                end
            end
        end
    end
    
    figure
    compare_trace_stack_grid(...
        spikes_grids,5,1,[],0,{'25','50','100'},1)

end
