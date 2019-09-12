function [trials] = design_batch_trials(neurons,design_params)
%%
trials(1) =struct;
i_trial =0;
for i= 1:length(neurons)
    this_location = neurons(i).location;
    jchosen=1;% for fixed design
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
            case 'off-nuclei'
                
                radius=design_params.candidate_grid_params.max_radius;
                n_point=design_params.nlocs_per_neuron-1;
                angles = (1:n_point)/n_point;
                tmp_loc=[[0 0 0]; [sin(2*pi*angles)' cos(2*pi*angles)' zeros(n_point,1)]];
                design_params.chosen_locations = ones(size(tmp_loc,1),1)*neurons(i).location+tmp_loc;
                chosen_flag=true;
            case 'linesearch'
                % Determine the chosen locations here: 
                design_params.grid_type='chosen';
                % There should only be two neurons
                one_loc=neurons(1).location;two_loc=neurons(2).location;
                numgrid=design_params.nlocs_per_neuron*2;
                gap=(two_loc-one_loc)/(numgrid-3);
                design_params.chosen_locations=ones(numgrid,1)*one_loc+...
                    [-1:(numgrid-2)]'*gap;
                              chosen_flag=true;
            case 'chosen'
                chose_flag=true;
        end
        
        if chosen_flag
            this_trial_location= design_params.chosen_locations(jchosen,:);
            n_rep=design_params.repeat_number;
            jchosen=jchosen+1;
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
