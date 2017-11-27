function [not_terminated] = check_all_learned(neighbourhoods)

learned_groups = {'disconnected', 'alive'};
learned_indicators= zeros(length(neighbourhoods),1);

for i_neighbourhood= 1:length(neighbourhoods)
    this_neighbourhood=neighbourhoods(i_neighbourhood);
    number_of_disconnected_cells= sum(get_group_inds(this_neighbourhood,learned_groups{1}));
    number_of_alive_cells= sum(get_group_inds(this_neighbourhood,learned_groups{2}));
    learned_indicators(i_neighbourhood)= (length(this_neighbourhood.neurons) ...
        == (number_of_alive_cells+number_of_disconnected_cells) );
%     fprintf('Nghbh %d: %d learned; ', i_neighbourhood, ((number_of_alive_cells+number_of_disconnected_cells)/length(this_neighbourhood.neurons)));

end
not_terminated = 1-prod(learned_indicators);
%  fprintf('\n Not terminated.\n');
end

%% For debugging
% remaining_cells = this_neighbourhood.neurons(get_group_inds(this_neighbourhood,'connected'));
% 
% 
% remaining_cells(1).PR_params(12)
% remaining_cells(1).truth
% 

