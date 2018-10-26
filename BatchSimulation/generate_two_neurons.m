function [neurons] = generate_two_neurons(simulation_params)
n=2;
neurons(n)=struct;
shift_tmp=[0 0 0];
for i= 1:n % since there are only two neurons
    neurons(i).truth=struct;
    neurons(i).truth.shift=shift_tmp;
    neurons(i).truth.optical_gain=0.01+(rand([1 1]))*0.04;
    neurons(i).truth.gain=neurons(i).truth.optical_gain;
    
    if simulation_params.batch.single_connection
        if i>1
            neurons(i).truth.PR=0;
        else
            neurons(i).truth.PR=0.9;
        end
    else
        neurons(i).truth.PR=0.2+unifrnd(0,1)*0.8;
    end
    if i>1
        switch simulation_params.batch.location
            case 'x'
                neurons(i).truth.location=[simulation_params.batch.location_radius 0 0];
            case 'y'
                neurons(i).truth.location=[0 simulation_params.batch.location_radius 0];
            case 'z'
                neurons(i).truth.location=[0 0 simulation_params.batch.location_radius];
            case 'a'
                a1=unifrnd(0,1)*2*pi;a2=unifrnd(0,1)*2*pi;
                neurons(i).truth.location=simulation_params.batch.location_radius*...
                    [sin(a1)*cos(a2) cos(a1)*cos(a2) sin(a2)];
        end
    else
        
        neurons(i).truth.location=[0 0 0];
    end
    %     one_neuron.truth.PR=gamma_truth(i_cell);
    if simulation_params.batch.delay_indicator
        neurons(i).truth.delay_mean=(rand(1)-0.5)*40+20;
        neurons(i).truth.delay_var=(rand(1)-0.5)*20+15;
    else
        neurons(i).truth.delay_mean=0;
        neurons(i).truth.delay_var=0;
    end
    
    neurons(i).truth.shape=[];
    neurons(i).fluorescence= []; % need to generate some fluorescence level
    neurons(i).cell_ID=i;
    
    neurons(i).location=neurons(i).truth.location+neurons(i).truth.shift;
    % Draw shapes:
end

