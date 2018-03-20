function plot_voltage_clamp_spatial(datafilename,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    batches = varargin{1};
else
    batches = [];
end

load(datafilename)

group_names = experiment_setup.group_names;
num_neighbourhoods = size(neighbourhoods,1);
num_batches = size(neighbourhoods,2);
if isempty(batches)
    batches = 2:num_batches;
end

if isfield(experiment_setup,'reproduced')
    experiment_queries = experiment_setup.records.queries;
end
% load batches and collect trials

    
for k = 1:num_neighbourhoods
    
    clear trials
    for i = batches
        
        experiment_query = experiment_queries(k,i);
        
        for j = 1:length(group_names)

            if isfield(experiment_query,group_names{j}) && ...
                    isfield(experiment_query.(group_names{j}),'trials') && ~isempty(experiment_query.(group_names{j}).trials)
                if ~exist('trials','var')
                    trials = experiment_query.(group_names{j}).trials;
                else
                    trials = [trials experiment_query.(group_names{j}).trials];
                end
            end

        end
        

    end
    
    for i = 1:length(neighbourhoods(k,end).neurons)
        if ~strcmp(neighbourhoods(k,end).neurons(i).group_ID{end},'secondary')
            this_ind = find([experiment_setup.neurons.cell_ID] == neighbourhoods(k,end).neurons(i).cell_ID);
            experiment_setup.neurons(this_ind).group_ID = neighbourhoods(k,end).neurons(i).group_ID{end};
        end
    end
    
    figure
    [vclamp_map, psc_time_map, color_map, linewidth_map, cell_map] = build_vclamp_grid(experiment_setup,trials,1);
    assignin('base','cell_map',cell_map)
    plot_timeseries_map(vclamp_map,Inf,1,0,color_map,linewidth_map,cell_map,psc_time_map);
    title(['Neighbourhood ' num2str(k)])
%     hold on
    

end
