function [mpp] = draw_samples(...
    trials_locations, trials_power, pi_target, background_rate,...
    cell_params,current_template,  funcs,    delay_params, stim_threshold,time_max)

   n_cell=length(cell_params);
    stimuli_size=zeros(size(trials_locations,1),n_cell);
   

for l = 1:size(trials_locations,1)
    for m = 1:size(trials_locations,2)
        if isnan(trials_locations(l,m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,trials_locations(l,m))*trials_power(l,m))';
        end
    end
end
n_trial =size(trials_locations,1);
mpp=struct();
mu_bg = 1/background_rate;

for i_trial = 1:n_trial
    mpp(i_trial).times=[];mpp(i_trial).assignments=[];
    mpp(i_trial).locations=trials_locations(i_trial,:);
    mpp(i_trial).power=trials_power(i_trial,:);
    
    for i_cell = 1:n_cell
        k=stimuli_size(i_trial,i_cell);
        if k > stim_threshold
        params_sim=cell_params(i_cell);
        
        stim = current_template*k;
        [~, spikes]  =lif_glm_sim(stim,params_sim,funcs);
%         plot(V_vect)
        if sum(spikes)>0
            if delay_params.type==2
                shape=(delay_params.mean^2)/(delay_params.std^2);
                scale = delay_params.mean/shape;
                delay_vec=round(gamrnd(shape,scale,[sum(spikes) 1]));
            else
                delay_vec=zeros([sum(spikes) 1]);
            end
            spikes_delay =find(spikes)+delay_vec';
            for i = 1:sum(spikes)
                if rand(1) < cell_params(i_cell).gamma
                % censoring at time max:
                if spikes_delay(i)<time_max
                    mpp(i_trial).times=[mpp(i_trial).times spikes_delay(i)];
                    mpp(i_trial).assignments=[mpp(i_trial).assignments i_cell];
                    
                end
                end
                
            end
        end
        end
    end
    % add background event:
    R = exprnd(mu_bg);
    while R < time_max
        mpp(i_trial).times=[mpp(i_trial).times max(1,round(R))];
        mpp(i_trial).assignments=[mpp(i_trial).assignments 0];
        R = R+exprnd(mu_bg);
    end
     
end
fprintf('%d trials', i_trial);
fprintf(' simulated;\n');

end

