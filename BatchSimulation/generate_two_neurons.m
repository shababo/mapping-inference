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
        neurons(i).truth.PR=rand(1)*0.4+0.6;
    else
        neurons(i).truth.PR= (rand(1)>0.5)*(rand(1)*0.4+0.6);
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


