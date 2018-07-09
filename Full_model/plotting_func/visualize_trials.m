function [figure_handle]=visualize_trials(figure_handle, experiment_query,neighbourhood,experiment_setup,varargin)
%%
if ~isempty(varargin) == 1
    show_change =varargin{1};
else
    show_change=false;
end
%%
neurons=experiment_setup.neurons;
% connected
color_set= {'k' 'g' 'w' 'r' 'y' 'b'};
alpha_threshold=0.1;
cell_IDs= [neurons(:).cell_ID];
size_scale =200;

subplot_dim= [0.1 0.1 0.55 0.83 ; 0.7 0.1 0.28 0.83];%l b w h


%% Visualize the neurons

neuron_coord=reshape([neurons(:).location],[3 length(neurons)])';
% color_set= {'r' 'b' 'g' 'w'};
itpr=zeros(length(neurons),2);
mean_pr=zeros(length(neurons),1);
current_groups = zeros(length(neurons),1);

for i_cell = 1:length(neurons)
    switch neurons(i_cell).group_ID{end}
        case 'undefined'
            col_id=1;
        case 'connected'
            col_id=2;
        case 'disconnected'
            col_id=4;
        case 'secondary'
            col_id = 5;
    end
    itpr(i_cell,1)=neurons(i_cell).posterior_stat(end).PR.lower_quantile;
    itpr(i_cell,2)=neurons(i_cell).posterior_stat(end).PR.upper_quantile;
    
    current_groups(i_cell)=col_id;
    mean_pr(i_cell)=neurons(i_cell).posterior_stat(end).PR.mean;
end
alpha_levels = (max(neuron_coord(:,3))-neuron_coord(:,3))/range(neuron_coord(:,3)) ;
alpha_levels(alpha_levels<alpha_threshold)=alpha_threshold;
%% Use the values in neighbourhood for the current update:
all_ID=[neurons(:).cell_ID];
if ~isempty(neighbourhood)
    neigh_IDs=[neighbourhood.neurons(:).cell_ID];
for i = 1:length(neighbourhood.neurons)
    cell_ID=neighbourhood.neurons(i).cell_ID;
    i_cell = find(all_ID==cell_ID);
    switch neighbourhood.neurons(i).group_ID{end}
        case 'undefined'
            col_id=1;
        case 'connected'
            col_id=2;
        case 'disconnected'
            col_id=4;
        case 'secondary' % do not change groupd id for secondary
            col_id = current_groups(i_cell);
    end
    itpr(i_cell,1)=neighbourhood.neurons(i).posterior_stat(end).PR.lower_quantile;
    itpr(i_cell,2)=neighbourhood.neurons(i).posterior_stat(end).PR.upper_quantile;
    
    current_groups(i_cell)=col_id;
    mean_pr(i_cell)=neighbourhood.neurons(i).posterior_stat(end).PR.mean;
    
end
end
% Obtain the values in the previous batch
if show_change
    pre_itpr=zeros(length(neurons),2);
    pre_mean_pr=zeros(length(neurons),1);
    pre_groups = zeros(length(neurons),1);
    
    for i = 1:length(neighbourhood.neurons)
        cell_ID=neighbourhood.neurons(i).cell_ID;
        i_cell = find(all_ID==cell_ID);
        switch neighbourhood.neurons(i).group_ID{end-1}
            case 'undefined'
                col_id=1;
            case 'connected'
                col_id=2;
            case 'disconnected'
                col_id=4;
            case 'secondary' % do not change groupd id for secondary
                col_id = current_groups(i_cell);
        end
        pre_itpr(i_cell,1)=neighbourhood.neurons(i).posterior_stat(end-1).PR.lower_quantile;
        pre_itpr(i_cell,2)=neighbourhood.neurons(i).posterior_stat(end-1).PR.upper_quantile;
        
        pre_groups(i_cell)=col_id;
        pre_mean_pr(i_cell)=neighbourhood.neurons(i).posterior_stat(end-1).PR.mean;
        
    end
    
end
%% Obtain trials:
group_names=fieldnames(experiment_setup.groups);
%
trials_by_group=struct;
for i_group =1:length(group_names)
    if isfield(experiment_query,group_names{i_group})
        these_trials=experiment_query.(group_names{i_group});
        if ~isempty(these_trials)
            loc_tmp=zeros(0,3);
            event_counts=[];
            for i_trial = 1:length( these_trials.trials)
                
                loc_tmp =[loc_tmp;these_trials.trials(i_trial).locations];
                if isfield(these_trials.trials(i_trial),'event_times')
                    event_counts =[event_counts; length( these_trials.trials(i_trial).event_times)];
                end
            end
            [loc_unique,ia,ic] = unique(loc_tmp,'rows');
            
            trials_by_group.(group_names{i_group})=struct;
            trials_by_group.(group_names{i_group}).locations = loc_unique;
            alpha_tmp=(max(neuron_coord(:,3))-loc_unique(:,3))/range(neuron_coord(:,3)) ;
            alpha_tmp(alpha_tmp<alpha_threshold)=alpha_threshold;
            alpha_tmp(alpha_tmp>1)=1;
            trials_by_group.(group_names{i_group}).alphas= alpha_tmp;
            event_flag= zeros(size(loc_unique,1),1);
            if isfield(these_trials.trials(i_trial),'event_times')
                for i_unique = 1:size(loc_unique,1)
                    event_flag(i_unique)=sum(event_counts(ic==i_unique));
                end
            end
            trials_by_group.(group_names{i_group}).event_flag=event_flag;
        end
    end
end

%%
% figure(figure_handle);
figure(1);

h(1)=subplot(1,2,1);
for i_cell = 1:length(neurons)
    this_color= color_set{current_groups(i_cell)};
    scatter(neuron_coord(i_cell,1),neuron_coord(i_cell,2),'SizeData', size_scale*mean_pr(i_cell),...
        'MarkerFaceColor',this_color,'MarkerEdgeColor',this_color,...
        'MarkerFaceAlpha', alpha_levels(i_cell),'MarkerEdgeAlpha', alpha_levels(i_cell))
    hold on;
    
    if strcmp(experiment_setup.experiment_type,'simulation')
        if (neurons(i_cell).truth.PR>0)
            scatter(neuron_coord(i_cell,1),neuron_coord(i_cell,2),'SizeData', size_scale*neurons(i_cell).truth.PR,...
                'MarkerEdgeColor',color_set{6},'LineWidth',2)%,...
            %'MarkerFaceAlpha',0)
            hold on;
        end
    end
    
    
end
% show boundary of neighbourhood:
% for i_neighbourhood = 1:length(neighbourhoods)
%     bound_tmp=neighbourhoods(i_neighbourhood).boundary;
% plot(bound_tmp(:,1),bound_tmp(1,2)*ones(2,1),'--',...
%     'Color','k','LineWidth',1)
% hold on;
% plot(bound_tmp(:,1),bound_tmp(2,2)*ones(2,1),'--',...
%     'Color','k','LineWidth',1)
% hold on;
% plot(bound_tmp(2,1)*ones(2,1),bound_tmp(:,2),'--',...
%     'Color','k','LineWidth',1)
% hold on;
% plot(bound_tmp(1,1)*ones(2,1),bound_tmp(:,2),'--',...
%     'Color','k','LineWidth',1)
% hold on;
% end

for i_group =1:length(group_names)
    if isfield(experiment_query,group_names{i_group})
        
        these_trials=experiment_query.(group_names{i_group});
        if ~isempty(these_trials)
            these_locations=trials_by_group.(group_names{i_group}).locations;
            these_alphas=trials_by_group.(group_names{i_group}).alphas;
            this_color=color_set{i_group};
            for i_loc = 1:size(these_locations,1)
                if trials_by_group.(group_names{i_group}).event_flag(i_loc)>0
                    marker='s';lwd=1;
                else
                    marker='x';lwd=3;
                end
                scatter(these_locations(i_loc,1),these_locations(i_loc,2),marker,...
                    'MarkerFaceColor',this_color,'MarkerEdgeColor',this_color,...
                    'MarkerFaceAlpha',these_alphas(i_loc),'MarkerEdgeAlpha',these_alphas(i_loc),...
                    'LineWidth',lwd)
                hold on;
            end
        end
    end
end
xlabel('x (um)');
ylabel('y (um)');
title('Cell map');


%%


h(2)= subplot(1,2,2);
[~, ps]=sort(itpr(:,2)+normrnd(0,1e-4,[length(neurons) 1]));
ranks= 1:length(neurons);
ranks(ps)=ranks;

for i_cell = 1:length(neurons)
    this_color= color_set{current_groups(i_cell)};
    plot(itpr(i_cell,:),[ranks(i_cell) ranks(i_cell)] ,...
        'Color',this_color,'LineWidth',1.5)
    hold on;
    scatter(mean_pr(i_cell),ranks(i_cell),...
        'MarkerFaceColor',this_color,'MarkerEdgeColor',this_color);
    hold on;
    
    if strcmp(experiment_setup.experiment_type,'simulation')
        if (neurons(i_cell).truth.PR>0)
            scatter(neurons(i_cell).truth.PR,ranks(i_cell),'Marker','^',...
                'MarkerFaceColor',color_set{6},'MarkerEdgeColor',color_set{5});
            hold on;
        end
    end
    
    if show_change
        cell_ID=neurons(i_cell).cell_ID;
        if sum(neigh_IDs==cell_ID)>0
            
            this_color= color_set{pre_groups(i_cell)};
            plot(pre_itpr(i_cell,:),[ranks(i_cell)-0.3 ranks(i_cell)-0.3] ,...
                'Color',this_color,'LineWidth',0.5)
            hold on;
            scatter(pre_mean_pr(i_cell),ranks(i_cell)-0.3,'SizeData',20,...
                'MarkerFaceColor',this_color,'MarkerEdgeColor',this_color);
            hold on;
        end
    end
end


plot(experiment_setup.groups.undefined.regroup_func_params.disconnected_threshold*ones(1,2),...
    [0 length(neurons)+1],'--',...
    'Color',color_set{4},'LineWidth',1.5)
hold on;
% plot(experiment_setup.groups.undefined.regroup_func_params.connected_threshold*ones(1,2),...
%     [0 length(neurons)+1],'--',...
%     'Color',color_set{2},'LineWidth',1.5)
% hold on;
ylim([0 length(neurons)+1]);
title('PR posterior');
set(h(1), 'position', subplot_dim(1,:) );
set(h(2), 'position', subplot_dim(2,:) );

%

% end

