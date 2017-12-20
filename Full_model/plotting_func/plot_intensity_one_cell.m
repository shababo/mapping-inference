function [figure_handle] = plot_intensity_one_cell(this_cell_ID,this_query,neighbourhood,experiment_setup,figure_handle,rescaled)

%%
target_cell_color='b';
overall_color='k';
event_color='r';
%% Extract all trials that are related to this cell


% current_group=neighbourhood.neurons([neighbourhood.neurons(:).cell_ID]==this_cell_ID).group_ID{end}

% find trials that stimulate this cell
% NOTE: we are missing cells that target its neighbours.
num_of_trials = length(this_query.trials);
related_trial_flag = zeros(num_of_trials,1);
for i_trial = 1:num_of_trials
    related_trial_flag(i_trial)=ismember(this_cell_ID,this_query.trials(i_trial).cell_IDs);
end


% related_trials = struct;
related_trial_index=find(related_trial_flag);
for i=1:sum(related_trial_flag)
    related_trials(i)=get_intensity_one_trial(this_query.trials(related_trial_index(i)),...
        neighbourhood,experiment_setup);
end

%% Plot:

% varargin for color and other graph specs
xgrid = (1:length(related_trials(1).estimated_intensity))/20;

% find the scale for intensity (in order to stack the intensity on y axis)
ymax=zeros(sum(related_trial_flag),1);
stim=zeros(sum(related_trial_flag),1); % for ordering the trials on the y axis
for i=1:sum(related_trial_flag)
    trial_intensity=related_trials(i);
    ymax(i)=max(trial_intensity.estimated_intensity);
    
    i_cell_this_trial=find([trial_intensity.stim_neurons(:).cell_ID]==this_cell_ID);
    stim(i)=trial_intensity.stim_neurons(i_cell_this_trial).stimulation;
end
[~,y_order]=sort(stim);

% DRAW ALL TRIALS (INTENSITY AND EVENTS)
if ~rescaled
    y_tick=cell([0 0]);
    for i = 1:length(y_order)
        y_tick{y_order(i)}= num2str(related_trials(i).trial_ID);
    end
    %
    
    figure(figure_handle)
    for i=1:sum(related_trial_flag)
        trial_intensity=related_trials(i);
        
        plot(xgrid,trial_intensity.estimated_intensity/max(ymax) + y_order(i)-1,'LineWidth',3,'Color',overall_color)
        hold on;
        if ~isempty(trial_intensity.event_times)
            for i_event = 1:length(trial_intensity.event_times)
                line(xgrid([trial_intensity.event_times(i_event),trial_intensity.event_times(i_event)]),[0.2,0.8]+y_order(i)-1,...
                    'LineStyle',':','LineWidth',3,'Color',event_color)
                hold on;
            end
        end
        i_cell_this_trial=find([trial_intensity.stim_neurons(:).cell_ID]==this_cell_ID);
        
        temp_intensity=trial_intensity.stim_neurons(i_cell_this_trial...
            ).intensity*trial_intensity.stim_neurons(i_cell_this_trial).PR;
        plot(xgrid,temp_intensity/max(ymax)+ y_order(i)-1,'LineWidth',3,'Color',target_cell_color)
        hold on;
    end
    %
    LegendHandels(1)= plot(nan,nan,'LineWidth',3,'Color',overall_color);
    hold on;
    legend_names{1}='Overall';
    LegendHandels(2)= plot(nan,nan,'LineStyle',':','LineWidth',3,'Color',event_color);
    hold on;
    legend_names{2}='Events';
    LegendHandels(3) =  ...
        plot(nan,nan,'LineWidth',3,'Color',target_cell_color);
    hold on;
    legend_names{3}=['Neuron', ' ', num2str(this_cell_ID)];
    
    legend(LegendHandels,legend_names,'Location','northeast');
    %
    xlabel('Time (ms)', 'FontSize',14)
    yticks(sort(y_order))
    yticklabels(y_tick)
    ylabel('Trial', 'FontSize',14)
    title(['Trials #', ' ',mat2str([related_trials.trial_ID])],'FontSize',14)
    
else % DRAW RESCALED EVENT TIMES
    %
    color_map = cool(length(y_order));
    all_events =[related_trials(:).scaled_times];
    
    figure(figure_handle)
    freq=histogram(all_events,'FaceColor','g');
    hold on;
    ymax=max(freq.Values);
    xmax=max(all_events)+1;
    for i=1:sum(related_trial_flag)
        trial_intensity=related_trials(i);
        
        if ~isempty(trial_intensity.event_times)
            for i_event = 1:length(trial_intensity.event_times)
                line([trial_intensity.scaled_times(i_event),trial_intensity.scaled_times(i_event)],[0.2,ymax/2],...
                    'LineStyle',':','LineWidth',3,'Color',color_map(y_order(i),:))
                hold on;
            end
        end
      
    end
    
    xlim([0 xmax]);
    ylim([0 ymax]);
     xlabel('Cumulative Intensity', 'FontSize',14)
    ylabel('Counts', 'FontSize',14)
    title(['Rescaled Events for neuron', ' ',num2str(trial_intensity.trial_ID)],'FontSize',14)
   
    %
end
