function [mpp]=extract_mpp(trials)
%trials=experiment_query_this_group.trials;

mpp=struct;
for i_trial = 1:length(trials)
   mpp(i_trial).event_times=trials(i_trial).event_times; 
    
end

end
