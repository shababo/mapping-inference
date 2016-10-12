function traces_by_location = view_hillel_data_by_location(file,power)

load(file,'sweeps','ExpStruct')

% num_trials = max(good_trials);
num_trials = length(sweeps);
stim_order = ExpStruct.MP285struct.stimOrder;
stim_power_order = ExpStruct.MP285struct.stimPowerOrder;
stim_locations = ExpStruct.MP285struct.StimPoints;
num_locations = size(stim_locations,1);
stim_order = stim_order(1:num_trials);
stim_power_order = stim_power_order(1:num_trials);

traces_by_location = cell(num_locations,2);

good_trials = zeros(1,num_trials);

for i = 1:num_trials
    if all(ExpStruct.stims{i}{5} == 0)
        good_trials(i) = 1;
    end
end


for i = 1:num_locations
    
%     good_trials_bool = 1:num_trials >= min(good_trials) & 1:num_trials <= max(good_trials);
    
    if isfield(ExpStruct.MP285struct,'stimPowerOrder')
        disp('Using strongest stim power...')
        if isnan(power)
            these_sweeps = find(good_trials & stim_order == i & stim_power_order ~= 4);
        else
            these_sweeps = find(good_trials & stim_order == i & stim_power_order == power);
        end
        sweeps_loc = sweeps(these_sweeps);
    else
        disp('Only one stim power...')
        sweeps_loc = sweeps(good_trials & stim_order == i);
    end
    traces_by_location{i,1} = zeros(length(sweeps_loc),size(sweeps{1}(2900:end,:),1));
    traces_by_location{i,2} = zeros(length(sweeps_loc),size(sweeps{1}(2900:end,:),1));
    
    for j = 1:length(sweeps_loc)
        
        traces_by_location{i,1}(j,:) = sweeps_loc{j}(2900:end,1);
        traces_by_location{i,2}(j,:) = sweeps_loc{j}(2900:end,2);
        
    end
end

try
    
    pia_normalized_pos = bsxfun(@minus,stim_locations,ExpStruct.MP285struct.piaPosition');
    [pia_vert_positions, pos_ordering] = sort(pia_normalized_pos(:,1));
    traces_by_location = traces_by_location(pos_ordering,:);
    vert_positions = pia_normalized_pos(pos_ordering,:);
    
    vert_positions
    
catch e
    disp('Error computing target positions')
    isfield(ExpStruct.MP285struct,'piaPosition')
    [vert_positions, pos_ordering] = sort(stim_locations(:,1));
    traces_by_location = traces_by_location(pos_ordering,:);
    vert_positions = stim_locations(pos_ordering,:);
end


try
    patched_cell_positions = bsxfun(@minus,[ExpStruct.MP285struct.cell1Position'; ExpStruct.MP285struct.cell2Position'], ExpStruct.MP285struct.piaPosition');

    patched_cell_positions
catch e
    disp('Error computing cell positions')
    isfield(ExpStruct.MP285struct,'cell1Position')
    isfield(ExpStruct.MP285struct,'cell2Position')
end

figure;
imagesc(squareform(pdist(vert_positions)))
    
figure
compare_trace_stack_grid({traces_by_location(:,1), traces_by_location(:,2)},Inf,1,[],0,{'raw','detected events'})
figure
compare_trace_stack_grid({traces_by_location(:,1), traces_by_location(:,2)},2,1,[],1,{'raw','detected events'})
figure
compare_trace_stack_grid_overlap({traces_by_location(:,1), traces_by_location(:,2)},Inf,1,[],0,{'L4','L5'},1)

