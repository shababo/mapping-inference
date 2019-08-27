function [] = visualize_fitted_trials(trials, neurons,plot_params)
ref_size=plot_params.ref_size;

if ~isfield(plot_params, 'by_neuron')
    plot_params.by_neuron=false;
end
if plot_params.prediction
    chosen_field = 'predicted';
else
      chosen_field = 'fitted';
end

figure(1)
tmp_loc=reshape([trials(:).locations],3,[])';
[locations_unique, ~, ic] = unique(tmp_loc,'rows');
for i_cell = 1:length(neurons)
    tmp=neurons(i_cell).location;
    scatter3(tmp(1),tmp(2),tmp(3),ref_size*2,  'MarkerFaceColor','red','MarkerEdgeColor','red');
    hold on;
    if plot_params.by_neuron
        text(tmp(1),tmp(2), tmp(3),num2str(i_cell));
    end
end

for i = 1:size(locations_unique,1)
    tmp=locations_unique(i,:);
    % Size scaled with power
    tmp_size = ref_size/2;
    scatter3(tmp(1),tmp(2), tmp(3), tmp_size,'MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerFaceAlpha',0.1);
    hold on;
    if ~plot_params.by_neuron
        text(tmp(1),tmp(2), tmp(3),num2str(i));
    end
end
hold off;


% Visualize the trials
if plot_params.by_neuron
    for i_neuron = 1:length(neurons)
        figure(i+1)
        these_trials =trials;
        plot_params.vertical_gap = 0.1;
        stim_size = zeros(length(these_trials),1);
        for i = 1:length(these_trials)
            stim_size(i)=these_trials(i).(chosen_field).stim(i_neuron+1);
        end
        [~, tmp]=sort(stim_size); % Ascending
        rankings = 1:length(these_trials);
        rankings(tmp)=rankings;
        visualize_fitted_trial_multiple(these_trials, rankings,plot_params)
        title(['Neuron ' num2str(i)])
    end
else
    for i = 1:size(locations_unique,1)
        figure(i+1)
        these_trials =trials(ic==i);
        plot_params.vertical_gap = 0.1;
        [~, tmp]=sort([these_trials(:).power_levels]); % Ascending
        rankings = 1:length(these_trials);
        rankings(tmp)=rankings;        
        visualize_fitted_trial_multiple(these_trials, rankings,plot_params)
        title(['Location ' num2str(i)])
    end
end



