% load 'work_spiking_data_0720.mat'
% load 'current_template.mat'

% load('spike_detect_revisit_about_to_fit_0810.mat','all_detection_grids')
load('data/current_template.mat')
load('data/all_detection_grids.mat')

%%
trial_length = 1500;
downres_rate = 20;
trial_length_downres = trial_length/downres_rate;

stims_t_norm = -1.0*norm_average_current(1:trial_length);
stims_t_norm_downres = downsample(stims_t_norm,downres_rate);
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

count = 1;

for i = 1:num_powers
    
    for l = 1:num_repeats
        for j = 1:sqrt(num_spatial_pos)
            for k = 1:sqrt(num_spatial_pos)
                stims_t_downres(count,:) = stims_t_norm_downres * powers(i);
                stims_x_vec(count,stims_x_value(count)) = 1;
                count = count + 1;
                
            end
        end
    end
end
%%
num_cells = length(all_detection_grids);

spikes_downres = cell(num_cells,1);
spatial_inits = zeros(11,11,num_cells);
all_detection_grids_downres = cell(size(all_detection_grids));

spike_count_locations = zeros(num_cells,sqrt(num_spatial_pos),sqrt(num_spatial_pos));

for m = 1:num_cells
    

    % stim_t = zeros(size(stims_x,1),1500);
    % % stim_t(100:199) = 1;
    % % stims_t = repmat(stim_t,size(stims_x,1),1);

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
                    
                    offset = 0;
                    spike_count_locations(m,j-offset,k-offset) = spike_count_locations(m,j-offset,k-offset)...
                        + sum(spikes_downres{m}(count,:));
                    count = count + 1;
                    
                end 
            end
        end
    end

end


%%
% 
cells_to_run = [1 2 3 4 5 6 7];
cells_to_run = 4;

g_vals = .5:.1:1.0;
% g_vals = .1
g_likelihoods = zeros([length(g_vals) length(cells_to_run)]);
num_params = 2 + num_spatial_pos;

fits = zeros(length(cells_to_run),length(g_vals),num_params);

voltages_power_grids = cell(3,length(cells_to_run),length(g_vals));
spikes_grids = cell(3,length(cells_to_run),length(g_vals));

spikes_per_location_data = cell(3,length(cells_to_run));
spikes_per_location_sim = cell(3,length(cells_to_run),length(g_vals));
first_spike_latency_data = cell(3,length(cells_to_run));
first_spike_latency_sim = cell(3,length(cells_to_run),length(g_vals));
spikes_per_location_var_data = cell(3,length(cells_to_run));
spikes_per_location_var_sim = cell(3,length(cells_to_run),length(g_vals));
first_spike_latency_var_data = cell(3,length(cells_to_run));
first_spike_latency_var_sim = cell(3,length(cells_to_run),length(g_vals));

for cell_choice_i = 1:length(cells_to_run);
    
    cell_choice = cells_to_run(cell_choice_i)
    
    spikes = spikes_downres{cell_choice}';
    % spikes = [spikes; zeros(25,num_trials)];

    count = 1;
    spike_locs = [];
    for j = 1:sqrt(num_spatial_pos)
        for k = 1:sqrt(num_spatial_pos)

            if spike_count_locations(cell_choice,j,k) > 0
                spike_locs = [spike_locs count];
            end
            count = count + 1;
        end
    end



    num_trials_good = num_trials - num_repeats*num_powers*(num_spatial_pos - length(spike_locs));
    bad_trials = [];
    for i = 1:num_trials
        if ~any(stims_x(i,1) == spike_locs)
            bad_trials = [bad_trials i];
        end
    end
    stims_t_downres_trunc = stims_t_downres;
    spikes_trunc = spikes;
    stims_t_downres_trunc(bad_trials,:) = [];
    spikes_trunc(:,bad_trials) = [];
    stims_x_trunc = stims_x;
    stims_x_trunc(bad_trials,:) = [];
    
    figure
    for m = 1:3

        spikes_per_location_data{m,cell_choice_i} = zeros(size(all_detection_grids_downres{cell_choice}{m}));
        first_spike_latency_data{m,cell_choice_i} = zeros(size(all_detection_grids_downres{cell_choice}{m}));
        spikes_per_location_var_data{m,cell_choice_i} = zeros(size(all_detection_grids_downres{cell_choice}{m}));
        first_spike_latency_var_data{m,cell_choice_i} = zeros(size(all_detection_grids_downres{cell_choice}{m}));


        for i = 1:size(all_detection_grids_downres{cell_choice}{m},1)
            for j = 1:size(all_detection_grids_downres{cell_choice}{m},2)

                spikes_per_location_data{m,cell_choice_i}(i,j) = mean(sum(all_detection_grids_downres{cell_choice}{m}{i,j}/70,2));
                spikes_per_location_var_data{m,cell_choice_i}(i,j) = var(sum(all_detection_grids_downres{cell_choice}{m}{i,j}/70,2));

                first_spikes_data = [];

                for k = 1:size(all_detection_grids_downres{cell_choice}{m}{i,j},1)

                    first_spikes_data = [first_spikes_data find(all_detection_grids_downres{cell_choice}{m}{i,j}(k,5:end),1,'first')];

                end

                first_spike_latency_data{m,cell_choice_i}(i,j) = mean(first_spikes_data);
                first_spike_latency_var_data{m,cell_choice_i}(i,j) = var(first_spikes_data);

            end
        end

        subplot(3,3,m)
        imagesc(spikes_per_location_data{m,cell_choice_i})
        colorbar
        title(['DATA: Cell ' num2str(cell_choice) ': Spikes/Location'])

        subplot(3,3,m + 3)
        imagesc(spikes_per_location_var_data{m,cell_choice_i})
        colorbar
        title(['DATA: Cell ' num2str(cell_choice) ': Spikes/Location Var'])

        subplot(3,3,m + 6)
%         imagesc(first_spike_latency_data{m,cell_choice_i})
        pcolor([first_spike_latency_data{m,cell_choice_i} nan(11,1); nan(1,11+1)]);
         shading flat;
         set(gca, 'ydir', 'reverse');
         caxis([0 60])
        colorbar
        title(['DATA: Cell ' num2str(cell_choice) ': First Spike Time Mean'])

%         subplot(4,3,m + 9)
% %         imagesc(first_spike_latency_var_data{m,cell_choice_i})
%         pcolor([first_spike_latency_var_data{m,cell_choice_i} nan(11,1); nan(1,11+1)]);
%          shading flat;
%          set(gca, 'ydir', 'reverse');
%         colorbar
%         title(['DATA: Cell ' num2str(cell_choice) ': First Spike Time Var'])

    end


    
    % max_ll = 

    for g_i = 1:length(g_vals)

        g = g_vals(g_i)
        
        [expg_hyperpol,expg_rheo,expg_stim]=gconv_multidim(stims_t_downres_trunc',spikes_trunc,g);

        full_stim_mat = zeros(trial_length_downres*num_trials_good,num_spatial_pos);

        for i = 1:num_trials_good
            full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x_trunc(i,1))...
                = expg_stim(:,i); 
        end

        stim_scale = 1/100;
        full_stim_mat = full_stim_mat*stim_scale;


        full_stim_mat_nozero = full_stim_mat(:,spike_locs);


        link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
        derlink = @(mu) exp(mu)./(exp(mu)-1);
        invlink = @(resp) log(1 + exp(resp));
        F = {link, derlink, invlink};

        [betahat_conv,~,stats_conv]=glmfit([expg_hyperpol(:) full_stim_mat_nozero],spikes_trunc(:),'binomial','link',F);

        disp(['Threshold: ' num2str(-1.0*betahat_conv(1))])
        disp(['Reset: ' num2str(betahat_conv(2))])

        spatial_filt_fit = zeros(num_spatial_pos,1);
        count = 1;
        for i = 1:num_spatial_pos
            if any(i == spike_locs)
                spatial_filt_fit(i) = betahat_conv(count+2);
                count = count + 1;
            end
        end
        
        fits(cell_choice_i,g_i,1:2) = betahat_conv(1:2);
        fits(cell_choice_i,g_i,3:end) = spatial_filt_fit;

        % compute ll
        [expg_hyperpol,expg_rheo,expg_stim]=gconv_multidim(stims_t_downres',spikes,g);
        full_stim_mat = zeros(trial_length_downres*num_trials,num_spatial_pos);

        for i = 1:num_trials_good
            full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x(i,1))...
                = expg_stim(:,i); 
        end

        stim_scale = 1/100;
        full_stim_mat = full_stim_mat*stim_scale;

        this_mu = betahat_conv(2)*expg_hyperpol(:) + betahat_conv(1) + full_stim_mat*spatial_filt_fit;
        this_lambda = log(exp(this_mu) + 1);
        g_likelihoods(g_i,cell_choice_i) = -sum(this_lambda) + sum(this_lambda.*spikes(:));

    

        figure; imagesc(reshape(spatial_filt_fit,sqrt(num_spatial_pos),sqrt(num_spatial_pos))')
        title(['Cell ' num2str(cell_choice) ', g = ' num2str(g) ': Spatial Filter, V_{th} = ' num2str(-1.0*betahat_conv(1)) ', V_{res} = ' num2str(betahat_conv(2))])

        %% output spikes


        %%%DEFINE PARAMETERS
        dt=1; %time step ms
        t_end=75; %total run time ms
        V_spike=70; %value to draw a spike to, when cell spikes

        E_L = 0;

        V_reset = betahat_conv(2);
        V_th = -betahat_conv(1);
        spatial_filt = zeros(num_spatial_pos,1);
        count = 1;
        for i = 1:num_spatial_pos
            if any(i == spike_locs)
                spatial_filt(i) = betahat_conv(count+2);
                count = count + 1;
            end
        end

        %%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
        t_vect=0:dt:t_end; 



        %%
        %INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m
        spikes_sim = zeros(t_end,num_trials);

        clear I_e_vect_mat
        
        voltages_grid = cell(11,11);
        spikes_grid = cell(11,11);

        for j = 1:num_trials %loop over different I_Stim values
            
            V_vect=zeros(1,length(t_vect));
            V_plot_vect=zeros(1,length(t_vect));
            V_plot_vect2=zeros(1,length(t_vect));
            lambda = zeros(1,length(t_vect));

            i=1; %index denoting which element of V is being assigned

            V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
            V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V
            V_plot_vect2(i)=V_vect(i);

            NumSpikes=0; %holds number of spikes that have occurred

            tao=exprnd(1);
            lambda(i)=log(exp(V_vect(i)-V_th)+1);
            last_spike = 1;

            for t=dt:dt:t_end %loop through values of t in steps of df ms        
                %V_inf = E_L + I_e_vect(i)*R_m;

                V_vect(i+1) = V_vect(i) + ...
                    (E_L-V_vect(i) + stims_t_downres(j,i)*spatial_filt(stims_x(j,1)))*g; %Euler's method
                lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);

                %if statement below says what to do if voltage crosses threshold
                if sum(lambda(last_spike+1:i+1))>tao %cell spiked
                    V_plot_vect2(i+1) = V_vect(i+1);
                    V_vect(i+1)=V_reset; %set voltage back to V_reset
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
        %             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
                else %voltage didn't cross threshold so cell does not spike
                    V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
                    V_plot_vect2(i+1) = V_vect(i+1);
        %             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
                end
                i=i+1; %add 1 to index,corresponding to moving forward 1 time step
            end    

            [ind1,ind2] = ind2sub([11 11],stims_x(j,1));
            voltages_grid{ind2,ind1} = [voltages_grid{ind2,ind1}; V_plot_vect];
            spikes_grid{ind2,ind1} = [spikes_grid{ind2,ind1}; spikes_sim(:,j)'];

        end

        %%


        for i = 1:3
            voltages_power_grids{i,cell_choice_i,g_i} = cell(11,11);
            for j = 1:11
                for k = 1:11
                    voltages_power_grids{i,cell_choice_i,g_i}{j,k} = voltages_grid{j,k}((i-1)*5+1:i*5,:);
                    spikes_grids{i,cell_choice_i,g_i}{j,k} = spikes_grid{j,k}((i-1)*5+1:i*5,:);
                end
            end
        end
        
        %%
        
        figure;compare_trace_stack_grid(voltages_power_grids(:,cell_choice_i,g_i),5,1,[],0,{'raw','detected events'})
        title(['Cell ' num2str(cell_choice) ', g = ' num2str(g) ': Output From LIF-GLM fit'])
        
        drawnow
        
        if g_i == 1
            figure;compare_trace_stack_grid(all_detection_grids_downres{cell_choice},5,1,[],0,{'raw','detected events'})
            title(['Cell ' num2str(cell_choice) ': Data (detected spikes)']);
        end
        
        %% summarize events

        figure
        for m = 1:3
            
            spikes_per_location_sim{m,cell_choice_i,g_i} = zeros(size(voltages_power_grids{m}));
            first_spike_latency_sim{m,cell_choice_i,g_i} = zeros(size(voltages_power_grids{m}));
            spikes_per_location_var_sim{m,cell_choice_i,g_i} = zeros(size(voltages_power_grids{m}));
            first_spike_latency_var_sim{m,cell_choice_i,g_i} = zeros(size(voltages_power_grids{m}));
            
            
            for i = 1:size(voltages_power_grids{m},1)
                for j = 1:size(voltages_power_grids{m},2)

                    spikes_per_location_sim{m,cell_choice_i,g_i}(i,j) = mean(sum(spikes_grids{m,cell_choice_i,g_i}{i,j},2));
                    spikes_per_location_var_sim{m,cell_choice_i,g_i}(i,j) = var(sum(spikes_grids{m,cell_choice_i,g_i}{i,j},2));

                    first_spikes_sim = [];
                    
                    for k = 1:size(spikes_grids{m}{i,j},1)

                        first_spikes_sim = [first_spikes_sim find(spikes_grids{m,cell_choice_i,g_i}{i,j}(k,5:end),1,'first')];

                    end
                    
                    first_spike_latency_sim{m,cell_choice_i,g_i}(i,j) = mean(first_spikes_sim);
                    first_spike_latency_var_sim{m,cell_choice_i,g_i}(i,j) = var(first_spikes_sim);
                    
                end
            end
            
            subplot(3,3,m)
            imagesc(spikes_per_location_sim{m,cell_choice_i,g_i})
            colorbar
            title(['Cell ' num2str(cell_choice) ', g = ' num2str(g) ': Spikes/Location'])
            
            subplot(3,3,m + 3)
            imagesc(spikes_per_location_var_sim{m,cell_choice_i,g_i})
            colorbar
            title(['Cell ' num2str(cell_choice) ', g = ' num2str(g) ': Spikes/Location Var'])
            
            subplot(3,3,m + 6)
%             imagesc(first_spike_latency_sim{m,cell_choice_i,g_i})
            pcolor([first_spike_latency_sim{m,cell_choice_i,g_i} nan(11,1); nan(1,11+1)]);
            shading flat;
            set(gca, 'ydir', 'reverse');
            caxis([0 60])
            colorbar
            title(['Cell ' num2str(cell_choice) ', g = ' num2str(g) ': First Spike Time Mean'])
            
%             subplot(4,3,m + 9)
% %             imagesc(first_spike_latency_var_sim{m,cell_choice_i,g_i})
%             pcolor([first_spike_latency_var_sim{m,cell_choice_i,g_i} nan(11,1); nan(1,11+1)]);
%             shading flat;
%             set(gca, 'ydir', 'reverse');
%             colorbar
%             title(['Cell ' num2str(cell_choice) ', g = ' num2str(g) ': First Spike Time Var'])
            
        end

    end
    
    figure; plot(g_vals,g_likelihoods(:,cell_choice_i)); title(['Cell ' num2str(cell_choice) ': g vs. LL'])
    
end
                
