load('data/current_template.mat')

%% build stimuli and set some params

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

%% set params for hard-thres lif for generating spikes

mu = [0 0];

Sigma = [1 0; 0 1];

x1 = -5:1:5; x2 = -5:1:5;
[X1,X2] = meshgrid(x1,x2);
spatial_filt_true = mvnpdf([X1(:) X2(:)],mu,Sigma);
spatial_filt_true = reshape(spatial_filt_true,length(x2),length(x1));

% put zeros in just to eff with it
spatial_filt_true(1,:) = 0;
spatial_filt_true(end,:) = 0;
spatial_filt_true(:,1) = 0;
spatial_filt_true(:,end) = 0;

V_th_true = 5; V_reset_true = -150; g_true = 0.05;
sf_gain = .5;
spatial_filt_true = spatial_filt_true*sf_gain;
spatial_filt_true = spatial_filt_true(:);


figure
imagesc(reshape(spatial_filt_true,sqrt(num_spatial_pos),sqrt(num_spatial_pos)));
colorbar
title('True Spatial Filt')

%%%DEFINE PARAMETERS
dt=1; %time step ms
t_end=75; %total run time ms
V_spike=70; %value to draw a spike to, when cell spikes

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 


%%
%INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m
spikes = zeros(t_end,num_trials);
voltages_grid = cell(11,11);
spikes_grid = cell(11,11);

for j = 1:num_trials %loop over different I_Stim values

    V_vect=zeros(1,length(t_vect));
    V_plot_vect=zeros(1,length(t_vect));
%     lambda = zeros(1,length(t_vect));

    i=1; %index denoting which element of V is being assigned

    V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V

    NumSpikes=0; %holds number of spikes that have occurred
    last_spike = 1;

    for t=dt:dt:t_end %loop through values of t in steps of df ms        
        %V_inf = E_L + I_e_vect(i)*R_m;

        V_vect(i+1) = V_vect(i) + ...
            -g_true*V_vect(i) + stims_t_downres(j,i)*spatial_filt_true(stims_x(j,1)); %Euler's method
%         lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        

        %if statement below says what to do if voltage crosses threshold
%         if sum(lambda(last_spike+1:i+1))>tao %cell spiked
        if V_vect(i+1) > V_th_true
            V_vect(i+1)=V_reset_true; %set voltage back to V_reset
            V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here

            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spikes(i,j)=1;
            if last_spike == i+1
                disp('two spikes in a row!!')
                return
            end
            last_spike = i+1;
%             tao=exprnd(1);
%             lambda(1:i)=0;
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        else %voltage didn't cross threshold so cell does not spike
            V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        end
        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
    end    

    [ind1,ind2] = ind2sub([11 11],stims_x(j,1));
    voltages_grid{ind2,ind1} = [voltages_grid{ind2,ind1}; V_plot_vect];
    spikes_grid{ind2,ind1} = [spikes_grid{ind2,ind1}; spikes(:,j)'];

end



spikes_grids = cell(3,1);
voltages_grids = cell(3,1);

for i = 1:3
    spikes_grids{i} = cell(11,11);
    voltages_grids{i} = cell(11,11);
    for j = 1:11
        for k = 1:11
            spikes_grids{i}{j,k} = spikes_grid{j,k}((i-1)*5+1:i*5,:)*70;
            voltages_grids{i}{j,k} = voltages_grid{j,k}((i-1)*5+1:i*5,:);
        end
    end
end

figure;compare_trace_stack_grid(spikes_grids,5,1,[],0,{'raw','detected events'})



spikes_per_location_true = cell(3,1);
first_spike_latency_true = cell(3,1);


figure
for m = 1:3

    spikes_per_location_true{m} = zeros(size(voltages_grids{m}));
    first_spike_latency_true{m} = zeros(size(voltages_grids{m}));



    for i = 1:size(voltages_grids{m},1)
        for j = 1:size(voltages_grids{m},2)

            spikes_per_location_true{m}(i,j) = mean(sum(spikes_grids{m}{i,j},2))/70;

            first_spikes_true = [];

            for k = 1:size(spikes_grids{m}{i,j},1)

                first_spikes_true = [first_spikes_true find(spikes_grids{m}{i,j}(k,:),1,'first')];

            end

            first_spike_latency_true{m}(i,j) = mean(first_spikes_true);

        end
    end

    subplot(2,3,m)
    imagesc(spikes_per_location_true{m})
    caxis([0 3])
    colorbar
    title(['Spikes/Location'])

    subplot(2,3,m + 3)
    pcolor([first_spike_latency_true{m} nan(11,1); nan(1,11+1)]);
    shading flat;
    set(gca, 'ydir', 'reverse');
    caxis([0 30])
    colorbar
    title(['First Spike Time Mean'])
end

colormap hot

%% fit this data

g_vals = [.01:.01:.1];
g_likelihoods = zeros([length(g_vals) 1]);
num_params = 2 + num_spatial_pos;
fits = zeros(length(g_vals),num_params);

spike_count_locations = zeros(sqrt(num_spatial_pos),sqrt(num_spatial_pos));

count = 1;
for i = 1:3
    for l = 1:5
        for j = 1:11
            for k = 1:11

                spike_count_locations(j,k) = spike_count_locations(j,k)...
                    + sum(spikes(:,count));
                count = count + 1;

            end 
        end
    end
end



count = 1;
spike_locs = [];
for j = 1:sqrt(num_spatial_pos)
    for k = 1:sqrt(num_spatial_pos)

        if spike_count_locations(j,k) > 0
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

stims_t_downres_noisey = stims_t_downres + normrnd(0,2.5,size(stims_t_downres));
% stims_t_downres_trunc = stims_t_downres_noisey;
stims_t_downres_trunc = stims_t_downres;
spikes_trunc = spikes;
stims_t_downres_trunc(bad_trials,:) = [];
spikes_trunc(:,bad_trials) = [];
stims_x_trunc = stims_x;
stims_x_trunc(bad_trials,:) = [];

% link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
% derlink = @(mu) exp(mu)./(exp(mu)-1);
% invlink = @(resp) log(1 + exp(resp));
link = @link_test;  %link = @(mu) mu + log(1-exp(-mu));
derlink = @derlink_test;
invlink = @invlink_test;

for g_i = 1:length(g_vals)

        g_tmp = g_vals(g_i)
        
        [expg_hyperpol,expg_rheo,expg_stim]=gconv_multidim(stims_t_downres_trunc',spikes_trunc,g_tmp);

        full_stim_mat = zeros(trial_length_downres*num_trials_good,num_spatial_pos);

        for i = 1:num_trials_good
            full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x_trunc(i,1))...
                = expg_stim(:,i); 
        end

        stim_scale = 1/1;
        full_stim_mat = full_stim_mat*stim_scale;


        full_stim_mat_trunc = full_stim_mat(:,spike_locs);

        F = {link, derlink, invlink};

        [betahat_conv,~,stats_conv]=glmfit([expg_hyperpol(:) full_stim_mat_trunc],spikes_trunc(:),'poisson','link',F);

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
%         spatial_filt_fit = betahat_conv(3:end);
        
        fits(g_i,1:2) = betahat_conv(1:2);
        fits(g_i,3:end) = spatial_filt_fit;

        % compute ll
        [expg_hyperpol,expg_rheo,expg_stim]=gconv_multidim(stims_t_downres_noisey',spikes,g_tmp);
        full_stim_mat = zeros(trial_length_downres*num_trials,num_spatial_pos);

        for i = 1:num_trials_good
            full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x(i,1))...
                = expg_stim(:,i); 
        end

        full_stim_mat = full_stim_mat*stim_scale;

        this_resp = betahat_conv(2)*expg_hyperpol(:) + betahat_conv(1) + full_stim_mat*spatial_filt_fit;
%         this_lambda = log(exp(this_resp) + 1);
        this_lambda = invlink_test(this_resp);
        g_likelihoods(g_i) = -sum(this_lambda) + sum(this_lambda.*spikes(:));

    

%         figure; imagesc(reshape(spatial_filt_fit,sqrt(num_spatial_pos),sqrt(num_spatial_pos))')
%         title(['g = ' num2str(g_tmp) ': Spatial Filter, V_{th} = ' num2str(-1.0*betahat_conv(1)) ', V_{res} = ' num2str(betahat_conv(2))])
        
end

figure; plot(g_vals,g_likelihoods)
%% re-output spikes

[~, g_mle_i] = max(g_likelihoods); 
g_test = g_vals(g_mle_i)
betahat_conv = fits(g_mle_i,:);

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
spatial_filt_sim = betahat_conv(3:end);

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 

% gain_sim = 100;
% V_reset_sim = V_reset*gain_sim;
% V_th_sim = V_th*gain_sim;
% spatial_filt_sim = spatial_filt*gain_sim;


%%
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
    lambda(i) = invlink_test(V_vect(i)-V_th_sim);
    last_spike = 1;

    for t=dt:dt:t_end %loop through values of t in steps of df ms        
        %V_inf = E_L + I_e_vect(i)*R_m;

        V_vect(i+1) = V_vect(i) + ...
            -g_test*V_vect(i) + stims_t_downres(j,i)*spatial_filt_sim(stims_x(j,1)); %Euler's method
%         lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
        lambda(i+1) = invlink_test(V_vect(i)-V_th_sim);

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

%%


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

