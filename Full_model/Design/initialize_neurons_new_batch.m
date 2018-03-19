function  [this_neighbourhood] =initialize_neurons_new_batch(this_neighbourhood)

% this_neighbourhood =neighbourhood;
number_of_cells=length(this_neighbourhood.neurons);
i_batch=this_neighbourhood.batch_ID;

if i_batch > 1
for i_cell = 1:number_of_cells
    this_neighbourhood.neurons(i_cell).PR_params(i_batch)=this_neighbourhood.neurons(i_cell).PR_params(i_batch-1);
    this_neighbourhood.neurons(i_cell).gain_params(i_batch)=this_neighbourhood.neurons(i_cell).gain_params(i_batch-1);
    this_neighbourhood.neurons(i_cell).delay_mu_params(i_batch)=this_neighbourhood.neurons(i_cell).delay_mu_params(i_batch-1);
    this_neighbourhood.neurons(i_cell).delay_sigma_params(i_batch)=this_neighbourhood.neurons(i_cell).delay_sigma_params(i_batch-1);
end


end

end
