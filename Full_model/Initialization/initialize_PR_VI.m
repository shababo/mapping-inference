function [variational_params]=initialize_PR_VI(variational_params,neurons,...
    trials,prior_info,inference_params,background_rate)
%%
%  reg_type='univariate';
reg_type='linear';
app_threshold =3;

n_cell=length(variational_params);
n_trial=length(trials);
initial_boundary_params = prior_info.prior_parameters.initial_boundary_params;
switch reg_type
    case 'univariate'
        for i_cell = 1:n_cell            
%             shifts=[variational_params(i_cell).shift_x.mean variational_params(i_cell).shift_y.mean...
%                 variational_params(i_cell).shift_z.mean];
%             this_adjusted_loc=neurons(i_cell).location-shifts;
            this_adjusted_loc=neurons(i_cell).location;
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
                PR_initial = max(0.01,PR_initial);
                variational_params(i_cell).PR.mean=log(PR_initial/(1-PR_initial));
            end
        end
    case 'linear'
        % create the outcome and design matrix:
        event_counts = zeros(length(trials),1);
        design_matrix = zeros(length(trials),n_cell);
        relevant_trials = zeros(n_cell,1);
        for i_trial = 1:length(trials)
            trials(i_trial).event_times=trials(i_trial).event_times( trials(i_trial).event_times < inference_params.event_range(2) &...
                trials(i_trial).event_times> inference_params.event_range(1));
            event_counts(i_trial) = length(trials(i_trial).event_times);
            for i_cell = 1:n_cell
                this_loc=neurons(i_cell).location;
                relevant_flag = false;
                for i_loc = 1:size(trials(i_trial).locations,1)
                    rel_position=trials(i_trial).locations(i_loc,:)-this_loc;
                    this_relevant_flag =check_in_boundary(rel_position,initial_boundary_params);
                    if  this_relevant_flag
                        relevant_flag =true;
                    end
                end
                if relevant_flag 
                    relevant_trials(i_cell)=true;
                end
                design_matrix(i_trial,i_cell)= relevant_flag;
            end
        end
        mdl = fitlm(design_matrix,event_counts-background_rate,'Intercept',false);
        betahat=mdl.Coefficients.Estimate;
        for i_cell = 1:n_cell
            if relevant_trials(i_cell)
                PR_initial=min(0.99,max(0.01,betahat(i_cell)));
                app_count= sum(design_matrix(:,i_cell));
                if app_count>app_threshold 
                variational_params(i_cell).PR.mean=log(PR_initial/(1-PR_initial));
                end
            end
            
        end
end