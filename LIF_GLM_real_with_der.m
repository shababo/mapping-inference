load('data/current_template_w_deriv.mat')

%%
rng(1234)

%% build stimuli and set some params

trial_length = 1500;
downres_rate = 20;
trial_length_downres = trial_length/downres_rate;

stims_t_norm = -1.0*norm_average_current(1:trial_length);
d_stims_t_norm = -1.0*d_norm_average_current(1:trial_length);
stims_t_norm_downres = downsample(stims_t_norm,downres_rate);
d_stims_t_norm_downres = downsample(d_stims_t_norm,downres_rate);

num_spatial_pos = 121;
powers = [25 50 100];
num_powers = length(powers);
num_repeats = 5;
num_trials = num_spatial_pos*num_repeats*num_powers;


stims_x = [repmat((1:num_spatial_pos)',num_repeats*num_powers,1) [25*ones(num_spatial_pos*num_repeats,1);...
            50*ones(num_spatial_pos*5,1); 100*ones(num_spatial_pos*num_repeats,1)]];

stims_x_value = repmat((1:num_spatial_pos)',num_repeats*num_powers,1);
stims_x_vec = zeros(num_trials,num_spatial_pos);


stims_t_downres = zeros(num_trials,trial_length_downres);
d_stims_t_downres = zeros(num_trials,trial_length_downres);

count = 1;

for i = 1:num_powers
    
    for l = 1:num_repeats
        for j = 1:sqrt(num_spatial_pos)
            for k = 1:sqrt(num_spatial_pos)
                stims_t_downres(count,:) = stims_t_norm_downres * powers(i);
                d_stims_t_downres(count,:) = d_stims_t_norm_downres * powers(i);
                stims_x_vec(count,stims_x_value(count)) = 1;
                count = count + 1;
                
            end
        end
    end
end

%%
load('data/all_detection_grids.mat')
%% 
num_cells = length(all_detection_grids);

spikes_downres = cell(num_cells,1);
spatial_inits = zeros(11,11,num_cells);
all_detection_grids_downres = cell(size(all_detection_grids));


for m = 1:num_cells
    

    spikes_downres{m} = zeros(num_trials,trial_length_downres);
    all_detection_grids_downres{m} = cell(1,3);
    count = 1;
    loc_count = 1;
    for i = 1:3
        all_detection_grids_downres{m}{i} = cell(11,11);
        for l = 1:5
            for j = 1:11
                for k = 1:11
                    if l == 1
                        all_detection_grids_downres{m}{i}{j,k} = zeros(5,75);
                    end
                    if l <= size(all_detection_grids{m}{i}{j,k},1)
                        spike_inds = ceil(find(all_detection_grids{m}{i}{j,k}(l,1:1500))/20);
                        spikes_downres{m}(count,spike_inds) = 1;
                        all_detection_grids_downres{m}{i}{j,k}(l,spike_inds) = 70;
                        if ~isempty(spike_inds)
                            spatial_inits(j,k,m) = spatial_inits(j,k,m) + 1;
                        end
                    else
                        spike_inds = [];
                    end
                    

                    count = count + 1;
                    
                end 
            end
        end
    end

end


spikes_per_location_data = cell(num_cells,1);
first_spike_latency_data = cell(num_cells,1);


for cell_i = 1:num_cells
    
    spikes_per_location_data{cell_i} = cell(3,1);
    first_spike_latency_data{cell_i} = cell(3,1);
    
    for m = 1:3
        spikes_per_location_data{cell_i}{m} = zeros(size(all_detection_grids_downres{cell_i}{m}));
        first_spike_latency_data{cell_i}{m} = zeros(size(all_detection_grids_downres{cell_i}{m}));



        for i = 1:size(all_detection_grids_downres{cell_i}{m},1)
            for j = 1:size(all_detection_grids_downres{cell_i}{m},2)

                spikes_per_location_data{cell_i}{m}(i,j) = mean(sum(all_detection_grids_downres{cell_i}{m}{i,j}/70,2));

                first_spikes_data = [];

                for k = 1:size(all_detection_grids_downres{cell_i}{m}{i,j},1)

                    first_spikes_data = [first_spikes_data find(all_detection_grids_downres{cell_i}{m}{i,j}(k,5:end),1,'first')];

                end

                first_spike_latency_data{cell_i}{m}(i,j) = mean(first_spikes_data);

            end
        end

    end
end

%% plot all stats

for i = 1:num_cells
    figure
    for m = 1:3
        
        subplot(2,3,m)
        imagesc(spikes_per_location_data{i}{m})
        caxis([0 3])
        colorbar
        title(['Spikes/Location'])

        subplot(2,3,m + 3)
        pcolor([first_spike_latency_data{i}{m} nan(11,1); nan(1,11+1)]);
        shading flat;
        set(gca, 'ydir', 'reverse');
        caxis([0 30])
        colorbar
        title(['First Spike Time Mean'])
    end
    colormap hot
end


%% fit this data

g_vals = [.01:.01:.1];
g_likelihoods = zeros([num_cells length(g_vals) 1]);
num_params = 2 + num_spatial_pos*2;
fits = zeros(num_cells,length(g_vals),num_params);

spike_locs = cell(num_cells,1);
for i = 1:num_cells
    
    spike_locs{i} = [];
    count = 1;
    for j = 1:11
        for k = 1:11

            if spikes_per_location_data{i}{1}(j,k) || ...
               spikes_per_location_data{i}{2}(j,k) || ...
               spikes_per_location_data{i}{3}(j,k)
                
                spike_locs{i} = [spike_locs{i} count];
                
            end
            
            count = count + 1;
                                       

        end 
    end    
end

% link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
% derlink = @(mu) exp(mu)./(exp(mu)-1);
% invlink = @(resp) log(1 + exp(resp));
link = @link_test;  %link = @(mu) mu + log(1-exp(-mu));
derlink = @derlink_test;
invlink = @invlink_test;

%%


for cell_i = 1:num_cells
    
    num_trials_good = num_trials - num_repeats*num_powers*(num_spatial_pos - length(spike_locs{cell_i}));
    bad_trials = [];
    for i = 1:num_trials
        if ~any(stims_x(i,1) == spike_locs{cell_i})
            bad_trials = [bad_trials i];
        end
    end

    % stims_t_downres = stims_t_downres + normrnd(0,2.5,size(stims_t_downres));
    stims_t_downres_trunc = stims_t_downres;
    d_stims_t_downres_trunc = d_stims_t_downres;
    spikes = spikes_downres{cell_i}';
    spikes_trunc = spikes;
    stims_t_downres_trunc(bad_trials,:) = [];
    d_stims_t_downres_trunc(bad_trials,:) = [];
    spikes_trunc(:,bad_trials) = [];
    stims_x_trunc = stims_x;
    stims_x_trunc(bad_trials,:) = [];


    for g_i = 1:length(g_vals)

            g_tmp = g_vals(g_i)

            [expg_hyperpol,expg_rheo,expg_stim,expg_dstim]=gconv_multidim(stims_t_downres_trunc',d_stims_t_downres_trunc',spikes_trunc,g_tmp);

            full_stim_mat = zeros(trial_length_downres*num_trials_good,num_spatial_pos);
            full_d_stim_mat = zeros(trial_length_downres*num_trials_good,num_spatial_pos);
            
            for i = 1:num_trials_good
                full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x_trunc(i,1))...
                    = expg_stim(:,i);
                full_d_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x_trunc(i,1))...
                    = expg_dstim(:,i);
            end

%             stim_scale = 1/1;
%             full_stim_mat = full_stim_mat*stim_scale;


            full_stim_mat_trunc = full_stim_mat(:,spike_locs{cell_i});
            full_d_stim_mat_trunc = full_d_stim_mat(:,spike_locs{cell_i});

            F = {link, derlink, invlink};

            [betahat_conv,~,stats_conv]=glmfit([expg_hyperpol(:) full_stim_mat_trunc full_d_stim_mat_trunc],spikes_trunc(:),'poisson','link',F);

            disp(['Threshold: ' num2str(-1.0*betahat_conv(1))])
            disp(['Reset: ' num2str(betahat_conv(2))])

            spatial_filt_fit = zeros(num_spatial_pos,1);
            d_spatial_filt_fit = zeros(num_spatial_pos,1);
            count = 1;
            for i = 1:num_spatial_pos
                if any(i == spike_locs{cell_i})
                    spatial_filt_fit(i) = betahat_conv(count+2);
                    count = count + 1;
                end
            end
            
            for i = 1:num_spatial_pos
                if any(i == spike_locs{cell_i})
                    d_spatial_filt_fit(i) = betahat_conv(count+2);
                    count = count + 1;
                end
            end
    %         spatial_filt_fit = betahat_conv(3:end);

            fits(cell_i,g_i,1:2) = betahat_conv(1:2);
            fits(cell_i,g_i,3:end) = [spatial_filt_fit; d_spatial_filt_fit];

            % compute ll
            [expg_hyperpol,expg_rheo,expg_stim,expg_dstim]=gconv_multidim(stims_t_downres',d_stims_t_downres',spikes,g_tmp);
            full_stim_mat = zeros(trial_length_downres*num_trials,num_spatial_pos);
            full_d_stim_mat = zeros(trial_length_downres*num_trials,num_spatial_pos);

            for i = 1:num_trials_good
                full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x(i,1))...
                    = expg_stim(:,i); 
                full_d_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x(i,1))...
                    = expg_dstim(:,i); 
            end

%             full_stim_mat = full_stim_mat*stim_scale;

            this_resp = betahat_conv(2)*expg_hyperpol(:) + betahat_conv(1) + full_stim_mat*spatial_filt_fit + full_d_stim_mat*d_spatial_filt_fit;
    %         this_lambda = log(exp(this_resp) + 1);
            this_lambda = invlink(this_resp);
            g_likelihoods(cell_i,g_i) = -sum(this_lambda) + sum(this_lambda.*spikes(:));



    %         figure; imagesc(reshape(spatial_filt_fit,sqrt(num_spatial_pos),sqrt(num_spatial_pos))')
    %         title(['g = ' num2str(g_tmp) ': Spatial Filter, V_{th} = ' num2str(-1.0*betahat_conv(1)) ', V_{res} = ' num2str(betahat_conv(2))])

    end

    figure; plot(g_vals,g_likelihoods(cell_i,:))
    title(sprintf('Cell %d',cell_i))

end
%% re-output spikes

all_sims_spikes_per_loc = cell(num_cells,1);
all_sims_first_spike_latency = cell(num_cells,1);
all_sims_spikes_grids = cell(num_cells,1);
all_sims_voltages_grids = cell(num_cells,1);

dt = 1;
t_end = 75;
t_vect=0:dt:t_end;
V_spike = 70;    

for cell_i = 1:num_cells
    
    [~, g_mle_i] = max(g_likelihoods(cell_i,:)); 
    g_test = g_vals(g_mle_i)
    betahat_conv = squeeze(fits(cell_i,g_mle_i,:));

    %%%DEFINE PARAMETERS
    V_reset_sim = betahat_conv(2);
    V_th_sim = -betahat_conv(1);
    % spatial_filt_sim = zeros(num_spatial_pos,1);
    % count = 1;
    % for i = 1:num_spatial_pos
    %     if any(i == spike_locs)
    %         spatial_filt_sim(i) = betahat_conv(count+2);
    %         count = count + 1;
    %     end
    % end
    spatial_filt_sim = betahat_conv(3:2+num_spatial_pos);
    d_spatial_filt_sim = betahat_conv(2+num_spatial_pos+1:end);




    
    %INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m
    spikes_sim = zeros(t_end,num_trials);

    voltages_grid_sim = cell(11,11);
    spikes_grid_sim = cell(11,11);

    for j = 1:num_trials %loop over different I_Stim values

        V_vect=zeros(1,length(t_vect));
        V_plot_vect=zeros(1,length(t_vect));
        lambda = zeros(1,length(t_vect));

        i=1; %index denoting which element of V is being assigned

        V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V

        NumSpikes=0; %holds number of spikes that have occurred

        tao=exprnd(1);
    %     lambda(i)=log(exp(V_vect(i)-V_th_sim)+1);
        lambda(i) = invlink(V_vect(i)-V_th_sim);
        last_spike = 1;

        for t=dt:dt:t_end %loop through values of t in steps of df ms        
            %V_inf = E_L + I_e_vect(i)*R_m;

            V_vect(i+1) = V_vect(i) + ...
                -g_test*V_vect(i) + stims_t_downres(j,i)*spatial_filt_sim(stims_x(j,1)) + d_stims_t_downres(j,i)*d_spatial_filt_sim(stims_x(j,1)); %Euler's method
    %         lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
            lambda(i+1) = invlink(V_vect(i)-V_th_sim);

            %if statement below says what to do if voltage crosses threshold
            if sum(lambda(last_spike+1:i+1))>tao %cell spiked
    %         if rand < lambda(i+1)
                V_vect(i+1)=V_reset_sim; %set voltage back to V_reset_sim
                V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here

                NumSpikes=NumSpikes+1; %add 1 to the total spike count
                spikes_sim(i,j)=1;
                if last_spike == i+1
                    disp('two spikes in a row!!')
                    return
                end
                last_spike = i+1;
                tao=exprnd(1);
    %             lambda(1:i)=0;
    %             lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
            else %voltage didn't cross threshold so cell does not spike
                V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
    %             lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
            end
            i=i+1; %add 1 to index,corresponding to moving forward 1 time step
        end    

        [ind1,ind2] = ind2sub([11 11],stims_x(j,1));
        voltages_grid_sim{ind2,ind1} = [voltages_grid_sim{ind2,ind1}; V_plot_vect];
        spikes_grid_sim{ind2,ind1} = [spikes_grid_sim{ind2,ind1}; spikes_sim(:,j)'];

    end

    


    spikes_grids_sim = cell(3,1);
    voltages_grids_sim = cell(3,1);

    for i = 1:3
        spikes_grids_sim{i} = cell(11,11);
        voltages_grids_sim{i} = cell(11,11);
        for j = 1:11
            for k = 1:11
                spikes_grids_sim{i}{j,k} = spikes_grid_sim{j,k}((i-1)*5+1:i*5,:)*70;
                voltages_grids_sim{i}{j,k} = voltages_grid_sim{j,k}((i-1)*5+1:i*5,:);
            end
        end
    end

    figure;compare_trace_stack_grid(spikes_grids_sim,5,1,[],0,{'raw','detected events'})

    spikes_per_location_sim = cell(3,1);
    first_spike_latency_sim = cell(3,1);


    figure
    for m = 1:3

        spikes_per_location_sim{m} = zeros(size(voltages_grids_sim{m}));
        first_spike_latency_sim{m} = zeros(size(voltages_grids_sim{m}));



        for i = 1:size(voltages_grids_sim{m},1)
            for j = 1:size(voltages_grids_sim{m},2)

                spikes_per_location_sim{m}(i,j) = mean(sum(spikes_grids_sim{m}{i,j},2))/70;

                first_spikes_sim = [];

                for k = 1:size(spikes_grids_sim{m}{i,j},1)

                    first_spikes_sim = [first_spikes_sim find(spikes_grids_sim{m}{i,j}(k,:),1,'first')];

                end

                first_spike_latency_sim{m}(i,j) = mean(first_spikes_sim);

            end
        end

        subplot(2,3,m)
        imagesc(spikes_per_location_sim{m})
        caxis([0 3])
        colorbar
        title(['Spikes/Location'])

        subplot(2,3,m + 3)
        pcolor([first_spike_latency_sim{m} nan(11,1); nan(1,11+1)]);
        shading flat;
        set(gca, 'ydir', 'reverse');
        caxis([0 30])
        colorbar
        title(['First Spike Time Mean'])
    end

    colormap hot
    
    all_sims_spikes_per_loc{cell_i} = spikes_per_location_sim;
    all_sims_first_spike_latency{cell_i} = first_spike_latency_sim;
    all_sims_spikes_grids{cell_i} = spikes_grids_sim;
    all_sims_voltages_grids{cell_i} = voltages_grids_sim;
    
    
    
end
