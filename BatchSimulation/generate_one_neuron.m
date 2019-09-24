function [one_neuron] = generate_one_neuron(simulation_params,prior_info)
prior_params=prior_info.prior_parameters;

one_neuron=struct;
one_neuron.truth=struct;
shift_tmp=[0 0 0];

tmp_params=prior_params;
tmp_params=rmfield(tmp_params,{'PR','background','shapes'});
one_neuron.truth.shift=shift_tmp;

one_neuron.truth.PR=1;
prior_sample= draw_samples_from_var_dist(tmp_params);

% one_neuron.truth.gain=prior_sample.gain;
one_neuron.truth.gain=prior_info.excite_inv( 0.04);

one_neuron.truth.shape=[];


one_neuron.truth.background=simulation_params.background_rate;
one_neuron.truth.location=[0 0 0];
%one_neuron.fluorescence= []; % need to generate some fluorescence level
one_neuron.location=one_neuron.truth.location+one_neuron.truth.shift;
one_neuron.cell_ID=1;

%     one_neuron.truth.PR=gamma_truth(i_cell);
if simulation_params.batch.delay_indicator
    one_neuron.truth.delay_mean=(rand(1)-0.5)*5+40;
    one_neuron.truth.delay_var=(rand(1)-0.5)*20+18;
else
    one_neuron.truth.delay_mean=0;
    one_neuron.truth.delay_var=0;
end
