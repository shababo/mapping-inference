function [neurons_both trials_both] = map_twoneurons_RD(simulation_params)

i_exp=5;
clear('trials_both', 'trials_c1', 'trials_c2');
n_trial = length(result_full_nrp(i_exp).spike_times_c1);
n_cell = 2;
trials_both(n_trial)=struct;
% trials_c1(n_trial)=struct;trials_c2(n_trial)=struct;

for i=1:n_trial 
    trials_both(i).power_levels=result_full_nrp(i_exp).spike_targ_power(i);
    trials_both(i).locations= result_full_nrp(i_exp).spike_targ_pos(i,:);
    
end
trials_c1=trials_both;trials_c2=trials_both;

for i=1:n_trial 
    etc1=result_full_nrp(i_exp).spike_times_c1(i);
    etc2=result_full_nrp(i_exp).spike_times_c2(i);
    if isnan(etc1)
       etc1=[];
    end
    if isnan(etc2)
       etc2=[];
    end
    
    trials_both(i).event_times = [etc1 etc2];
    trials_both(i).truth.assignments = [ones(length(etc1), 1)*1 ones(length(etc2), 1)*2];
    
    trials_c1(i).event_times = etc1;
    trials_c1(i).truth.assignments =ones(length(etc1), 1)*1;
    
    trials_c2(i).event_times = etc2;
    trials_c2(i).truth.assignments =ones(length(etc2), 1)*1;
end

clear('neurons_both')
neurons_both(2)=struct;
neurons_both(1).location=result_full_nrp(i_exp).c1_pos;
neurons_both(1).cell_ID=1;
neurons_both(2).location=result_full_nrp(i_exp).c2_pos;
neurons_both(2).cell_ID=2;
