function plot_voltage_clamp_spatial(datafilename)

load(datafilename)

group_names = experiment_setup.group_names;
num_neighbourhoods = size(neighbourhoods,1);
num_batches = size(neighbourhoods,2);

% load batches and collect trials

for k = 1:num_neighbourhoods
    
    clear trials
    for i = 2:num_batches
        
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
    
    figure
    [vclamp_map, psc_time_map, color_map] = build_vclamp_grid(experiment_setup,trials,1);
    plot_timeseries_map(vclamp_map,Inf,1,0,color_map,[],[],psc_time_map);
    title(['Neighbourhood ' num2str(k)])

end
