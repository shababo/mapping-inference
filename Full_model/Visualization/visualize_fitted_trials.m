function [] = visualize_fitted_trials(trials, neurons,plot_params)
%%
ref_size=plot_params.ref_size;

if ~isfield(plot_params, 'by_neuron')
    plot_params.by_neuron=false;
end
if plot_params.prediction
    chosen_field = 'predicted';
else
      chosen_field = 'fitted';
end
if ~isfield(plot_params, 'stim_scale')
     plot_params.stim_scale = @(x) x;
end
fig_num = randi(10000);
i = 0;
% figure(fig_num)
tmp_loc=reshape([trials(:).locations],3,[])';
[locations_unique, ~, ic] = unique(tmp_loc,'rows');
% for i_cell = 1:length(neurons)
%     tmp=neurons(i_cell).location;
%     scatter3(tmp(1),tmp(2),tmp(3),ref_size*2,  'MarkerFaceColor','red','MarkerEdgeColor','red');
%     hold on;
%     if plot_params.by_neuron
%         text(tmp(1),tmp(2), tmp(3),num2str(i_cell));
%     end
% end
% 
% for i = 1:size(locations_unique,1)
%     tmp=locations_unique(i,:);
%     % Size scaled with power
%     tmp_size = ref_size/2;
%     scatter3(tmp(1),tmp(2), tmp(3), tmp_size,'MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerFaceAlpha',0.1);
%     hold on;
%     if ~plot_params.by_neuron
%         text(tmp(1),tmp(2), tmp(3),num2str(i));
%     end
% end
% hold off;

% Find the maximum intensity and use it to standardize all intensities 
plot_params.max_intensity=1e-2;
for i_trial = 1:length(trials)
    this_trial =trials(i_trial);
    if plot_params.spike
        fits = this_trial.(chosen_field).intensity.spike;
    else
        fits =this_trial.(chosen_field).intensity.event;
    end
    plot_params.max_intensity = max(plot_params.max_intensity, max(max(fits)));
end

% Visualize the trials
if plot_params.by_neuron
    for i_neuron = 1:length(neurons)
       fig= figure(i_neuron+1+fig_num);
       plot_params.loc_indices=ic;
       
       if   plot_params.event_only % only show locations with recorded events 
       these_trials =trials(plot_params.evt_flag==1);
       plot_params.these_indices=plot_params.loc_indices(plot_params.evt_flag==1);
      
       else
       these_trials =trials;
       plot_params.these_indices=plot_params.loc_indices;
       end 
        plot_params.gap = 0.1;
        plot_params.colors=lines(size(locations_unique,1)); 
         stim_size = zeros(length(these_trials),1);
        
        for i = 1:length(these_trials)
            stim_size(i)=these_trials(i).(chosen_field).stim(i_neuron+1);
        end
        stim_size = log(stim_size);
%         [~, tmp]=unique(stim_size); % Ascending
%         
%         rankings = 1:length(these_trials);
%         rankings(tmp)=rankings;

        [covered_flags]=visualize_fitted_trial_multiple(these_trials, stim_size,plot_params);
        title(plot_params.main_title,'FontSize', plot_params.lab_size)
        yyaxis left
ylim([0 plot_params.time_max]);

yyaxis right
ylabel('Expected # events','FontSize', plot_params.lab_size)
ylim([0 plot_params.time_max]);
        if isfield(plot_params,'save_path')
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
            save_path = [plot_params.save_path plot_params.typename '_All.png'];
            saveas(fig,save_path)
        end
    end
else
%     for i = 1:size(locations_unique,1)
% these_trials =trials(ic==i);
%         fig=figure(i+1+fig_num);
% %         plot_params.vertical_gap = 0.1;
% % %         [~, tmp]=sort([these_trials(:).power_levels]); % Ascending
% % %         rankings = 1:length(these_trials);
% % %         rankings(tmp)=rankings;        
% 
%         plot_params.these_indices=plot_params.loc_indices(ic==i);
%         visualize_fitted_trial_multiple(these_trials,[these_trials(:).power_levels],plot_params)
%         title(['Location ' num2str(i)],'FontSize', plot_params.lab_size)
%         if isfield(plot_params,'save_path')
%             save_path = [plot_params.save_path '_' num2str(i+1) '.png'];
%             saveas(fig,save_path)
%         end
%     end
end

% set(gca, 'XScale',  plot_params.stim_scale)
% xlim([-3 5])


% if isfield(plot_params,'off_loc')
%     fig= figure(i_neuron+4+fig_num)
%        these_trials =trials;
%        plot_params.vertical_gap = 0.1;
%        dist_list = zeros(length(these_trials),1);
%         for i = 1:length(these_trials)
%             
%             rel_pos=these_trials(i).locations-neurons(1).location;
%             dist_list(i)= min( sqrt(sum(rel_pos(1:2).^2)),abs(rel_pos(3)));
%         end
%        visualize_fitted_trial_multiple(these_trials,dist_list,plot_params)
%       ylim([0,2.5])
%      ylabel('Trials ordered by deviation from xy-plane and z-axis');
% end
% 
%% Draw the coverage of the true spikes in groundtruth data:
if isfield(trials(1),'truth') 
    
    rec=plot_params.spike;
plot_params.spike =true;
plot_params.max_intensity=1e-2;
for i_trial = 1:length(trials)
    this_trial =trials(i_trial);
    fits = this_trial.(chosen_field).intensity.spike;
      plot_params.max_intensity = max(plot_params.max_intensity, max(max(fits)));
end

% Visualize the trials
    for i_neuron = 1:length(neurons)
       fig= figure(i_neuron+1+2*fig_num);
       plot_params.loc_indices=ic;
       
       if   plot_params.event_only % only show locations with recorded events
           these_trials =trials(plot_params.evt_flag==1);
           plot_params.these_indices=plot_params.loc_indices(plot_params.evt_flag==1);
       else
           these_trials =trials;
           plot_params.these_indices=plot_params.loc_indices;
       end
        plot_params.gap = 0.1;
        plot_params.colors=lines(size(locations_unique,1)); 
        
        spike_times = zeros(length(these_trials),1);
        
        for i = 1:length(these_trials)
            these_trials(i).event_times=these_trials(i).truth.spike_times;
           spike_times(i)=these_trials(i).truth.spike_times;
        end
%         spike_times_loc=zeros(max(ic),1);
%         for i=1:max(ic)
%             spike_times_loc(i)=mean(spike_times(ic==i));
%         end
%        spike_times_mean = spike_times_loc(ic);

        visualize_fitted_trial_multiple(these_trials, spike_times,plot_params);
        title(plot_params.main_title,'FontSize', plot_params.lab_size)
        
xlabel('Spike time')
yyaxis left
ylim([0 plot_params.time_max]);

yyaxis right
ylabel('Expected # events','FontSize', plot_params.lab_size)
ylim([0 plot_params.time_max]);


        if isfield(plot_params,'save_path')
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
            save_path = [plot_params.save_path plot_params.typename '_Spike.png'];
            saveas(fig,save_path)
        end
    end
    plot_params.spike=rec;
end 
%% Draw a spatial map: 
fig=figure(10+fig_num);

for i=1:length(neurons)
    this_loc = neurons(i).location;
scatter3(this_loc(1),this_loc(2),this_loc(3),plot_params.markerSize*3,"s",'filled',...
        'MarkerFaceColor', 'k',  'MarkerEdgeColor','k')
hold on;
end
% Map plot_params.evt_flag to locations:
evt_flag_loc = zeros(size(locations_unique,1),1);
for i= 1:size(locations_unique,1)
   evt_flag_loc(i)=max(plot_params.evt_flag(ic==i));
end

for i = 1:size(locations_unique,1)
    if evt_flag_loc(i)==1
    scatter3(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3),plot_params.markerSize,'filled',...
        'MarkerFaceColor', plot_params.colors(i,:),  'MarkerEdgeColor',plot_params.colors(i,:))
    hold on;
    if ~covered_flags(i) 
       text(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3), ['Loc ' num2str(i)]);
    end
    else
        scatter3(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3),plot_params.markerSize,...
         'MarkerEdgeColor',plot_params.colors(i,:))
    hold on;
    end
end
if isfield(plot_params,'save_path')
    save_path = [plot_params.save_path plot_params.typename '_Map.png'];
    saveas(fig,save_path)
end
%% Draw the trials at the problematic locations: 
if  plot_params.event_only
plot_params.fit_type='full_intensity';
for i = 1:size(locations_unique,1)
    plot_params.by_neuron=false;
    if ~covered_flags(i)
        these_trials =trials(ic==i);
        fig=figure(10+i+fig_num);
        plot_params.these_indices=plot_params.loc_indices(ic==i);
        %         plot_params.vertical_gap = 0.1;
        % %         [~, tmp]=sort([these_trials(:).power_levels]); % Ascending
        % %         rankings = 1:length(these_trials);
        % %         rankings(tmp)=rankings;
        visualize_fitted_trial_multiple(these_trials,[these_trials(:).power_levels],plot_params);
        title(['Location ' num2str(i)],'FontSize', plot_params.lab_size)
        
        ylim([0 5+plot_params.time_max]);
            if isfield(plot_params,'save_path')
                        mkdir([plot_params.save_path  plot_params.typename ]);
                save_path = [plot_params.save_path  plot_params.typename  '/Loc' num2str(i) '.png'];
                saveas(fig,save_path)
            end
    end
end
end
