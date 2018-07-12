 function [variational_params]=initialize_PR_VI(variational_params,neurons,trials,prior_info)
   
 
n_cell=length(variational_params);
n_trial=length(trials);
boundary_params = prior_info.prior_parameters.boundary_params;

for i_cell = 1:n_cell
    
    shifts=[variational_params(i_cell).shift_x.mean variational_params(i_cell).shift_y.mean...
        variational_params(i_cell).shift_z.mean];
    this_adjusted_loc=neurons(i_cell).location-shifts;
    relevant_trials=[];
    responsive_counts = 0;
    for i_trial = 1:length(trials)
        relevant_flag = false;
        for i_loc = 1:size(trials(i_trial).locations,1)
            rel_position=trials(i_trial).locations(i_loc,:)-this_adjusted_loc;
            this_relevant_flag =check_in_boundary(rel_position,boundary_params);
            if  this_relevant_flag
                relevant_flag =true;
            end
        end
        if relevant_flag
            relevant_trials = [relevant_trials i_trial];
            responsive_counts =  responsive_counts+(1-isempty(trials(i_trial).event_times));
        end
    end
    if ~isempty(relevant_trials)
        PR_initial=responsive_counts/length(relevant_trials);
        PR_initial = max(0.1,PR_initial);
    variational_params(i_cell).PR.mean=log(PR_initial/(1-PR_initial));
    end
end
