function [trials] = design_batch_trials(neurons,design_params)
%%
trials(1) =struct;
i_trial =0;
for i= 1:length(neurons)
    this_location = neurons(i).location;
    for j= 1:design_params.nlocs_per_neuron
        switch design_params.grid_type
            case 'random_2d'
                radius=design_params.candidate_grid_params.max_radius;
                this_radius = unifrnd(0,radius);
                this_angle = unifrnd(0,1);
            case 'random_3d'
                radius=design_params.candidate_grid_params.max_radius;
                this_radius = unifrnd(0,radius);
                this_angle = unifrnd(0,1);
                
                this_angle2=unifrnd(0,1);
            case 'ring'
                radius=design_params.candidate_grid_params.max_radius;
                n_ring=length(design_params.candidate_grid_params.number);
                i_ring = randsample(1:n_ring,1,true,design_params.candidate_grid_params.number);
                n_ring_grid=design_params.candidate_grid_params.number(i_ring);
                this_radius= (i_ring-1)*(radius/(n_ring-1));
                this_angle = randsample( 1:n_ring_grid,1)/n_ring_grid;
            case 'nuclei'
                this_radius=0;this_angle=0;
        end
        if strcmp(design_params.grid_type,'chosen')
            this_trial_location= design_params.chosen_locations(j,:);
            n_rep=design_params.repeat_number;
        else
            if j==1 & design_params.always_nucleus
                this_radius=0; n_rep=design_params.repeat_number_nucleus;
            else
                n_rep=design_params.repeat_number;
            end
            if strcmp(design_params.grid_type,'random_3d')
            this_trial_location=this_location+...
                this_radius*[sin(2*pi*this_angle)*sin(2*pi*this_angle2) cos(2*pi*this_angle)*sin(2*pi*this_angle2) cos(2*pi*this_angle2)];
            
            else
            this_trial_location=this_location+...
                this_radius*[sin(2*pi*this_angle) cos(2*pi*this_angle) 0];
                
            end
        end
        
        for i_pow = 1:length(design_params.power_levels)
            for i_rep=1:n_rep
                    i_trial = i_trial +1;
        
                   trials(i_trial).locations=this_trial_location;
                trials(i_trial).power_levels=design_params.power_levels(i_pow);
            end
        end
        
    end
end
