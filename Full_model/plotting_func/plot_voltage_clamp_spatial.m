function plot_voltage_clamp_spatial(experiment_setup,neighbourhood_ID,batch_IDs,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    
    dirpath = varargin{1};
 
else
    
    dirpath = experiment_setup.analysis_root;
end
linespecs = cell(0);    
group_names = experiment_setup.group_names;

% load batches and collect trials
for i = 1:length(batch_IDs)
    
    batchsavepath = [dirpath experiment_setup.exp_id ...
                    '_n' num2str(neighbourhood_ID)...
                    '_b' num2str(batch_IDs(i)) '_complete.mat'];
                
    load(batchsavepath)
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

[vclamp_map, psc_time_map, color_map] = build_vclamp_grid(experiment_setup,trials,3);
    

plot_timeseries_map(vclamp_map,Inf,1,0,color_map,[],[],psc_time_map);
