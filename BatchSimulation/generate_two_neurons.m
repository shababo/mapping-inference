function [neurons] = generate_two_neurons(simulation_params,prior_params)
n=2;
neurons(n)=struct;
shift_tmp=[0 0 0];

tmp_params=prior_params;
tmp_params=rmfield(tmp_params,{'PR','background','shapes'});

for i= 1:n
    
    
    prior_sample= draw_samples_from_var_dist(tmp_params);
    neurons(i).truth=struct;
    neurons(i).truth.shift=shift_tmp;
    neurons(i).truth.gain=prior_sample.gain;
    
    if i==1
        neurons(i).truth.PR=rand(1)*0.5+0.3;
    else
        neurons(i).truth.PR= (rand(1)>0.5)*(rand(1)*0.5+0.3);
    end
    
    %     one_neuron.truth.PR=gamma_truth(i_cell);
    if simulation_params.batch.delay_indicator
        neurons(i).truth.delay_mean=prior_sample.delay_mean;
        neurons(i).truth.delay_var=prior_sample.delay_var;
    else
        neurons(i).truth.delay_mean=0;
        neurons(i).truth.delay_var=0;
    end
    
    neurons(i).truth.shape=[];
    neurons(i).fluorescence= []; % need to generate some fluorescence level
    neurons(i).cell_ID=i;
    neurons(i).truth.location=[0 0 0];
    neurons(i).location=neurons(i).truth.location+neurons(i).truth.shift;
    % Draw shapes:
end



switch simulation_params.batch.location
    case 'x'
        neurons(2).truth.location=[simulation_params.batch.location_radius 0 0];
    case 'y'
        neurons(2).truth.location=[0 simulation_params.batch.location_radius 0];
    case 'z'
        neurons(2).truth.location=[0 0 simulation_params.batch.location_radius];
    case 'a'
        a1=unifrnd(0,1)*2*pi;a2=unifrnd(0,1)*2*pi;
        neurons(2).truth.location=simulation_params.batch.location_radius*...
            [sin(a1)*cos(a2) cos(a1)*cos(a2) sin(a2)];
end
  neurons(2).location=neurons(2).truth.location+neurons(2).truth.shift;
