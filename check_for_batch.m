function [batch_found,neighbourhood_ID] = check_for_batch(experiment_setup,neighbourhoods)

% check for batch on these neighborhoods
batch_found = 0;
for i = 1:length(neighbourhoods)
    neighbourhood = neighbourhoods(i);
    neighbourhood_ID = neighbourhood.neighbourhood_ID;
    batch_ID = neighbourhood.batch_ID + 1; % look for next batch
    check_file = [experiment_setup.phase_mask_dir experiment_setup.exp_id ...
                '_n' num2str(neighbourhood_ID)...
                '_b' num2str(batch_ID) '_batch_ready.mat'];
        
    batch_found = fileattrib(check_file);
    if batch_found
        break
    end
end

           
    