function [dat] = pre_processing_curves(single_patch_path,varargin)

load(single_patch_path); % should contain result_current,result_spikes

x_current=[];
spike_curve_cell_ids = [];
y_spike_mean=[];
y_spike_sd=[];  

for i_cell = 1:length(result_current)
    if ~isempty(result_current(i_cell).peak_current_means)
%     i_cell

    %     if i_cell ~= 11
%     if length(result_current(i_cell).these_powers) ==  length(result_spikes(i_cell).these_powers)
%         if min(result_current(i_cell).these_powers ==  result_spikes(i_cell).these_powers)==1 
            % NOTE: using these two criteria will leave only 4 cells.
            % It seems that not all the cells have matching power levels
            % between the two structures

        [~,i_c,i_s] = intersect(result_current(i_cell).these_powers,result_spikes(i_cell).these_powers);
            % make sure the power levels line with those in result_current and
            % result_spikes
            x_current =[ x_current result_current(i_cell).peak_current_means(i_c)];
            spike_curve_cell_ids = [spike_curve_cell_ids i_cell*ones(size(result_current(i_cell).peak_current_means(i_c)))];
            y_spike_mean=[y_spike_mean result_spikes(i_cell).spike_time_means(i_s)];
            y_spike_sd=[y_spike_sd result_spikes(i_cell).spike_time_jitter(i_s)];

%         end
%     end
    %     end
    end
end

dat.x_current=x_current;
dat.spike_curve_cell_ids = spike_curve_cell_ids ;
dat.y_spike_mean=y_spike_mean;
dat.y_spike_sd=y_spike_sd;
