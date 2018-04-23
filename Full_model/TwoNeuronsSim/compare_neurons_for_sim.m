function [results]=compare_neurons_for_sim(fitted_neurons,neurons,spike_curves)
results(length(neurons))=struct;
% 1) whether it learns the true parameters:
% 2) whether it predicts the event time (& spike time) well 
for i =1:length(neurons)
    results(i).delay_mean= fitted_neurons(i).delay_mean-neurons(i).delay_mean;
    results(i).delay_var= fitted_neurons(i).delay_var-neurons(i).delay_var;
    results(i).PR= fitted_neurons(i).PR-neurons(i).PR;
    results(i).gain= fitted_neurons(i).gain-neurons(i).gain;
    results(i).mpp(length(fitted_neurons(i).mpp))=struct;
    for j= 1:length(fitted_neurons(i).mpp)
        event_index=find(fitted_neurons(i).mpp(j).assignments==i);
        if isempty(event_index)
            results(i).mpp(j).event_times = [];
            results(i).mpp(j).spike_times = [];
            results(i).mpp(j).stimulation = [];
            
            results(i).mpp(j).pred_spike_times = [];
            results(i).mpp(j).pred_event_times = [];
            results(i).mpp(j).location=[];
            results(i).mpp(j).spike_dev = [];
            results(i).mpp(j).event_dev =[];
            
            
        else
            results(i).mpp(j).event_times = fitted_neurons(i).mpp(j).event_times(event_index);
            results(i).mpp(j).spike_times = fitted_neurons(i).mpp(j).spike_times(event_index);
            results(i).mpp(j).stimulation = fitted_neurons(i).mpp(j).stimulation(event_index);
            stim=results(i).mpp(j).stimulation*fitted_neurons(i).gain;
            results(i).mpp(j).pred_spike_times = predict_spike_time(stim,spike_curves);
            results(i).mpp(j).pred_event_times = results(i).mpp(j).pred_spike_times+fitted_neurons(i).delay_mean;
            results(i).mpp(j).location=fitted_neurons(i).mpp(j).location;
            results(i).mpp(j).spike_dev = results(i).mpp(j).pred_spike_times-results(i).mpp(j).spike_times;
            results(i).mpp(j).event_dev = results(i).mpp(j).pred_event_times-results(i).mpp(j).event_times;
            
        end
    end
end
