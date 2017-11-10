function [neurons] = generate_neurons(simulation_params)
%
        
% There should be options to read neuron locations 
%% Simulation some more neurons:


cell_locations=rand(simulation_params.number_of_cells,size(simulation_params.dimension,2));
for i_dimension= 1:size(simulation_params.dimension,2)
    cell_locations(:,i_dimension)=cell_locations(:,i_dimension)*range(simulation_params.dimension(:,i_dimension))+...
        simulation_params.dimension(1,i_dimension);
end

extra_locations = 2*pi*rand(simulation_params.siblings.number,size(simulation_params.dimension,2)); % only need the first two columns 
for i= 1:simulation_params.siblings.number 
    temp_cell=randsample(1:simulation_params.number_of_cells,1);
    extra_locations(i,:)=cell_locations(temp_cell,:)+ ...
        simulation_params.siblings.distance*[sin(extra_locations(i,1))*sin(extra_locations(i,2)),...
        cos(extra_locations(i,1))*sin(extra_locations(i,2)),cos(extra_locations(i,2))];
end

cell_locations=[cell_locations;extra_locations];

number_of_cells = size(cell_locations,1);


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

 \item neurons(cell\_ID)
    fluorescence
        \item location
        \item gain\_prior \note{calculated from fluorescence level and location}
        \item PR\_prior \note{the following are from \struct{prior\_info}} 
\item        V\_reset
\item        V\_thresh
\item        membrane\_resistance 
\item truth \note{for simulation}
        \end{itemize}
    
neurons=struct([]);
for i_cell = 1:number_of_cells 
    neurons(i_cell).fluorescence= []; % need to generate some fluorescence level 
    neurons(i_cell).V_reset= -1e4;
    neurons(i_cell).V_thresh=15;
    neurons(i_cell).membrane_resistance=0.02;
    neurons(i_cell).location=cell_locations(i_cell,:);
    neurons(i_cell).optical_gain=[];
    neurons(i_cell).PR=[];
    neurons(i_cell).truth=struct;
        neurons(i_cell).truth.V_reset= -1e4;
        neurons(i_cell).truth.V_thresh=15;
        neurons(i_cell).truth.membrane_resistance=0.02;
        neurons(i_cell).truth.location=cell_locations(i_cell,:);
        neurons(i_cell).truth.optical_gain=gain_truth(i_cell);
        neurons(i_cell).truth.PR=gamma_truth(i_cell);
end

