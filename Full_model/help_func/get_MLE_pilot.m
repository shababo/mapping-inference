function [gain_this_cell,gain_sample_mat,loglklh]=get_MLE_pilot(result_spikes_this_cell, spike_curves,varargin)
% result_spikes_this_cell contains info only about one cell
% call:
% for this_cell = 1:length(result_spikes)
%   [gain_mle(this_cell)]=get_MLE_pilot(result_spikes(this_cell), spike_curves);
% end

if ~isempty(varargin) && ~isempty(varargin{1})
    specs= varargin{1};
else
    specs=struct;
    specs.gain_gap=0.001;
    specs.low=0.001;
    specs.up=0.4;
    
    specs.Tmax=500;
    specs.power_ref=0;
    specs.background_rate=1e-6;% just to prevent singularity 
    specs.lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
    
end
% 
% %%
% result_spikes_this_cell=result_spikes(this_cell);
%   specs=struct;
%     specs.gain_gap=0.001;
%     specs.low=0.005;
%     specs.up=0.3;
%     
%     specs.Tmax=500;
%     specs.power_ref=4.79;
%     specs.background_rate=3e-4;
%     specs.lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;

%%
gain_sample_mat = specs.low:specs.gain_gap:specs.up; % grid of gains

number_of_trials = length(result_spikes_this_cell.power{1});
stim_size = result_spikes_this_cell.power{1}-specs.power_ref;
mpp=struct;
for i_trial = 1:number_of_trials
    mpp(i_trial).event_times = result_spikes_this_cell.spike_times{1}(i_trial);
    if isnan(mpp(i_trial).event_times)
        mpp(i_trial).event_times=[];
    elseif mpp(i_trial).event_times>specs.Tmax
        mpp(i_trial).event_times=specs.Tmax;
    end
end
S=size(gain_sample_mat,2);
result_spikes_this_cell=length(mpp);
% n_grid=size(prob_trace_full,2);

t_grid= 1:specs.Tmax;
loglklh=zeros(1,S);

for s=1:length(gain_sample_mat)
    gain_sample=gain_sample_mat(s);
    loglklh_vec = zeros(result_spikes_this_cell,1);
    for  i_trial = 1:result_spikes_this_cell
        
        stim_temp =stim_size(i_trial,:)';
        effective_stim= stim_temp.*gain_sample;
        
        % find the index of the actual stimulation in the vector:
        [~, stim_index]=min(abs(effective_stim - spike_curves.current));
        prob_this_trial=zeros(1,specs.Tmax);
        prob_this_trial(:)=...
            normpdf(t_grid,spike_curves.mean(stim_index),spike_curves.sd(stim_index));
        prob_this_trial=[specs.background_rate*ones(1, specs.Tmax); prob_this_trial];
        
        [loglklh_vec(i_trial)]=  specs.lklh_func(mpp(i_trial),prob_this_trial);
        
    end
    
    loglklh(s)=sum(loglklh_vec);
    
end
[~, i_gain]=max(loglklh);
gain_this_cell=gain_sample_mat(i_gain);

%%
% scatter(gain_sample_mat,loglklh)