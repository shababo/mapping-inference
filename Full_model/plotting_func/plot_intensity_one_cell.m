function [figure_handle] = plot_intensity_one_cell(this_cell_ID,this_query,neighbourhood,...
    experiment_setup,figure_handle,rescaled,varargin)

if ~isempty(varargin)
   plot_type='collapsed';
   stim_gap=varargin{1};
else
   plot_type='stacked';
end
if rescaled
   plot_type='rescaled'; 
end
%%
target_cell_color='b';
overall_color='k';
event_color='r';
%% Extract all trials that are related to this cell


% current_group=neighbourhood.neurons([neighbourhood.neurons(:).cell_ID]==this_cell_ID).group_ID{end};
current_group=this_query.trials(1).group_ID;
% find trials that stimulate this cell
% NOTE: we are missing cells that target its neighbours.
stim_size = get_stim_size(current_group,this_query.trials,neighbourhood);
i_this_cell = find([neighbourhood.neurons(:).cell_ID]==this_cell_ID);
effective_stim_this_cell=stim_size(:,i_this_cell)*neighbourhood.neurons(i_this_cell).gain_params(end).mean;
minimum_stim_threshold=experiment_setup.prior_info.induced_intensity.minimum_stim_threshold;

related_trial_flag = effective_stim_this_cell>(minimum_stim_threshold/3);

num_of_trials = length(this_query.trials);


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
switch plot_type
    case 'stacked'
        %%
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
                    'LineStyle','-','LineWidth',3,'Color',event_color)
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
    LegendHandels(2)= plot(nan,nan,'LineStyle','-','LineWidth',3,'Color',event_color);
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
    
    case 'rescaled'
    %%
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
    
    case 'collapsed'
    %%
    stim_bounds = cell([1 1]);
    stim_bounds{1}=[min(stim) min(stim)+stim_gap];
    i=2;
    while stim_bounds{i-1}(2) < max(stim)
        stim_bounds{i}(1)=min(stim(stim> stim_bounds{i-1}(2)));
        stim_bounds{i}(2)=stim_bounds{i}(1)+stim_gap;
        i=i+1;
    end
    
     y_tick=cell([0 0]);
     group_trials = cell([0 0]);
    for i = 1:length(stim_bounds)
        
        group_trials{i}= find(  stim>(stim_bounds{i}(1)-0.1) & stim<(stim_bounds{i}(2)+0.1));
        y_tick{i}= num2str( round(mean(stim(group_trials{i}))) );
    end
    %
    
    figure(figure_handle)
    for i=1:length(stim_bounds)
        for j = 1:length(group_trials{i})
            
            trial_intensity=related_trials(group_trials{i}(j));
            
            plot(xgrid,trial_intensity.estimated_intensity/max(ymax) +i-1,'LineWidth',1,'Color',overall_color)
            hold on;
            if ~isempty(trial_intensity.event_times)
                for i_event = 1:length(trial_intensity.event_times)
                    line(xgrid([trial_intensity.event_times(i_event),trial_intensity.event_times(i_event)]),[0.2,0.8]+i-1,...
                        'LineStyle','-','LineWidth',2,'Color',event_color)
                    hold on;
                end
            end
            i_cell_this_trial=find([trial_intensity.stim_neurons(:).cell_ID]==this_cell_ID);
            
            temp_intensity=trial_intensity.stim_neurons(i_cell_this_trial...
                ).intensity*trial_intensity.stim_neurons(i_cell_this_trial).PR;
            plot(xgrid,temp_intensity/max(ymax)+ i-1,'LineWidth',1,'Color',target_cell_color)
            hold on;
        end
        
    end
    %
    LegendHandels(1)= plot(nan,nan,'LineWidth',3,'Color',overall_color);
    hold on;
    legend_names{1}='Overall';
    LegendHandels(2)= plot(nan,nan,'LineStyle','-','LineWidth',3,'Color',event_color);
    hold on;
    legend_names{2}='Events';
    LegendHandels(3) =  ...
        plot(nan,nan,'LineWidth',3,'Color',target_cell_color);
    hold on;
    legend_names{3}=['Neuron', ' ', num2str(this_cell_ID)];
    
    legend(LegendHandels,legend_names,'Location','northeast');
    %
    xlabel('Time (ms)', 'FontSize',14)
    yticks(0:(length(stim_bounds)-1))
    yticklabels(y_tick)
    ylabel('Stimulation', 'FontSize',14)
    title(['Trials #', ' ',mat2str([related_trials.trial_ID])],'FontSize',14)
   
end
