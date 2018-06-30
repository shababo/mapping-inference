 function [variational_params, prior_params]=initialize_params_VI(neurons)
    n_cell=length(neurons);
    clear('variational_params')
    %variational_params(n_cell)=struct;
    for i_cell = 1:n_cell
        variational_params(i_cell)=neurons(i_cell).params(end);
    end
    prior_params=variational_params;
   