function [trials] = design_batch_trials(neurons,design_params)
% Design trials to stimulate each neuron in a given number of trials 
%   - stimulation location can be random or on designed grids near the
%   nuclei


% design_params.ntrials_per_neuron = 10;
% design_params.grid_type = 'ring';
% design_params.candidate_grid_params.number=[1 6 12];
% design_params.candidate_grid_params.max_radius=10;
% design_params.power_levels = 10:10:40;

% [neighbourhood.neurons(i_cell_nhood).posterior_stat(end).shift_x.mean ...
% neighbourhood.neurons(i_cell_nhood).posterior_stat(end).shift_y.mean ...
% neighbourhood.neurons(i_cell_nhood).posterior_stat(end).shift_z.mean];
trials(1) =struct;
radius=design_params.candidate_grid_params.max_radius;
i_trial =0;
for i= 1:length(neurons)
    this_location = neurons(i).location;
    for j= 1:design_params.nlocs_per_neuron
        i_trial = i_trial +1;
        trials(i_trial).power_levels=randsample(design_params.power_levels,1,true);
        
        switch design_params.grid_type
            case 'random'
                this_radius = unifrnd(0,radius);
                this_angle = unifrnd(0,1);
            case 'ring'
                n_ring=length(design_params.candidate_grid_params.number);
                i_ring = randsample(1:n_ring,1,true,design_params.candidate_grid_params.number);
                n_ring_grid=design_params.candidate_grid_params.number(i_ring);
                
                this_radius= (i_ring-1)*(radius/(n_ring-1));
                this_angle = randsample( 1:n_ring_grid,1)/n_ring_grid;
        end
        if j==1 & design_params.always_nucleus
            this_radius=0;
        end
        trials(i_trial).locations=this_location+...
            this_radius*[sin(2*pi*this_angle) cos(2*pi*this_angle) 0];
        if design_params.repeat_number > 1
            if j==1 & design_params.always_nucleus
            n_rep=design_params.repeat_number_nucleus;
            else
            n_rep=design_params.repeat_number;
            end
            for i=2:n_rep
                i_trial = i_trial +1;
                trials(i_trial).locations=this_location+...
            this_radius*[sin(2*pi*this_angle) cos(2*pi*this_angle) 0];
                 trials(i_trial).power_levels=randsample(design_params.power_levels,1,true);
            end
        end 
    end
end
