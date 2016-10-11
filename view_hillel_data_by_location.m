function traces_by_location = view_hillel_data_by_location(file)

load(file,'sweeps','ExpStruct')

num_trials = length(sweeps);
stim_order = ExpStruct.MP285struct.stimOrder;
stim_locations = ExpStruct.MP285struct.StimPoints;
num_locations = size(stim_locations,1);
stim_order = stim_order(1:num_trials);

traces_by_location = cell(num_locations,2);

for i = 1:num_locations
    
    sweeps_loc = sweeps(stim_order == i);
    traces_by_location{i,1} = zeros(length(sweeps_loc),size(sweeps{1}(2900:end,:),1));
    traces_by_location{i,2} = zeros(length(sweeps_loc),size(sweeps{1}(2900:end,:),1));
    
    for j = 1:length(sweeps_loc)
        
        traces_by_location{i,1}(j,:) = sweeps_loc{j}(2900:end,1);
        traces_by_location{i,2}(j,:) = sweeps_loc{j}(2900:end,2);
        
    end
end

figure
compare_trace_stack_grid({traces_by_location(:,1), traces_by_location(:,2)},10,1,[],0,{'raw','detected events'})
figure
compare_trace_stack_grid({traces_by_location(:,1), traces_by_location(:,2)},10,1,[],1,{'raw','detected events'})
figure
compare_trace_stack_grid_overlap({traces_by_location(:,1), traces_by_location(:,2)},10,1,[],0,{'L4','L5'},1)

