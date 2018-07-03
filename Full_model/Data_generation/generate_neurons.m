function [neurons] = generate_neurons(experiment_setup)
% To-dos:
% 1) simulate shifts
% 2) unify the names of fields

simulation_params = experiment_setup.sim;
% There should be options to read neuron locations
%% Simulation some more neurons:

if simulation_params.use_real_map
    if isempty(simulation_params.real_map_name)
        warning('Real cell maps not specified.')
        
    else
        if strcmp(experiment_setup.experiment_type,'simulation')
            cell_map_path =[experiment_setup.exp_root simulation_params.real_map_name];
        else
            cell_map_path =[experiment_setup.analysis_root simulation_params.real_map_name];
        end
        if ~exist(cell_map_path, 'file')
            warning(['Specified cell maps not found in ' cell_map_path '.'])
        else
            
            load(cell_map_path,'cell_locs')
            cell_locations=cell_locs;
            % Overwrite the dimension info:
            for i_dimension = 1:size(cell_locations,2)
                simulation_params.dimension(:,i_dimension)= [min(cell_locations(:,i_dimension)) max(cell_locations(:,i_dimension))];
                
            end
        end
    end
else
    cell_locations=rand(simulation_params.number_of_cells,size(simulation_params.dimension,2));
    for i_dimension= 1:size(simulation_params.dimension,2)
        cell_locations(:,i_dimension)=cell_locations(:,i_dimension)*range(simulation_params.dimension(:,i_dimension))+...
            simulation_params.dimension(1,i_dimension);
    end
end
number_of_cells = size(cell_locations,1);

extra_locations = 2*pi*rand(simulation_params.siblings.number,size(simulation_params.dimension,2)); % only need the first two columns
for i= 1:simulation_params.siblings.number
    temp_cell=randsample(1:number_of_cells,1);
    extra_locations(i,:)=cell_locations(temp_cell,:)+ ...
        simulation_params.siblings.distance*[sin(extra_locations(i,1))*sin(extra_locations(i,2)),...
        cos(extra_locations(i,1))*sin(extra_locations(i,2)),cos(extra_locations(i,2))];
end

cell_locations=[cell_locations;extra_locations];

funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
switch simulation_params.cell_params.type
    case 'Normal' % normal
        gamma_truth = (rand([number_of_cells 1])<simulation_params.connection_params.proportion).*(0.5+0.3*rand([number_of_cells 1]));
        gain_truth=0.02+(rand([number_of_cells 1])-0.5)*0.01;
    case 'Extreme gain' % normal gamma, extreme gains
        gain_bound_truth.up=0.03;gain_bound_truth.low=0.005;
        
        gamma_truth = (rand([number_of_cells 1])<simulation_params.connection_params.proportion).*(0.5+0.3*rand([number_of_cells 1]));
        n_connected=sum(gamma_truth>0);
        lognormal_temp=exp(normrnd(0,sqrt(1),[number_of_cells  1]));
        gain_truth= (gain_bound_truth.up-gain_bound_truth.low)*exp(-lognormal_temp)./(1+exp(-lognormal_temp))+...
            gain_bound_truth.low;
        gain_truth=max(gain_bound_truth.low+0.0005,gain_truth);
    case 'Weak gamma' % weak gamma, normal gains
        weak_gamma_proportion=0.5; % some gammas are small
        gamma_truth = (rand([number_of_cells  1])<simulation_params.connection_params.proportion).*(0.5+0.3*rand([number_of_cells  1]));
        n_connected=sum(gamma_truth>0);
        weak_index=randsample(find(gamma_truth>0), floor(weak_gamma_proportion*n_connected));
        gamma_truth(weak_index)=(0.15+0.1*rand([length(weak_index) 1]));
        gain_truth=0.02+(rand([number_of_cells  1])-0.5)*0.01;
    otherwise
        
        % throw a warning
end

%%
prior_params =experiment_setup.prior_info.prior_parameters;
%% Store all simulated parameters in the truth field, including the shape mean & kernel
%
neurons=struct([]);
for i_cell = 1:number_of_cells
    
    neurons(i_cell).truth=struct;
    %     neurons(i_cell).truth.V_reset= -1e4;
    %     neurons(i_cell).truth.V_thresh=15;
    %     neurons(i_cell).truth.membrane_resistance=0.02;
    neurons(i_cell).truth.location=cell_locations(i_cell,:);
    shift_tmp=[normrnd(prior_params.shift_x.mean,exp(prior_params.shift_x.log_sigma)),...
        normrnd(prior_params.shift_y.mean,exp(prior_params.shift_y.log_sigma)),...
        normrnd(prior_params.shift_z.mean,exp(prior_params.shift_z.log_sigma))];
    neurons(i_cell).truth.shift=shift_tmp;
    
    neurons(i_cell).truth.optical_gain=gain_truth(i_cell);
    neurons(i_cell).truth.PR=gamma_truth(i_cell);
    neurons(i_cell).truth.delay_mean=(rand(1)-0.5)*40+40;
    neurons(i_cell).truth.delay_var=(rand(1)-0.5)*20+15;
    tmp= draw_3D_GP(simulation_params.mesh_grid,1,prior_params.GP_params);
    neurons(i_cell).truth.shape= tmp.samples; % n_loc by 1 vector
           
    neurons(i_cell).fluorescence= []; % need to generate some fluorescence level
    %     neurons(i_cell).V_reset= -1e4;
    %     neurons(i_cell).V_thresh=15;
    %     neurons(i_cell).membrane_resistance=0.02;
    
    neurons(i_cell).cell_ID=i_cell;
    
    neurons(i_cell).location=neurons(i_cell).truth.location+neurons(i_cell).truth.shift;
    neurons(i_cell).cell_ID = i_cell;
    
end

