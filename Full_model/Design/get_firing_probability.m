function [prob_trace] = get_firing_probability(...
    linkfunc,...
    current_template,stim_unique,cell_params,delay_params)
%----------- Parameters of the template cell

num_stim_unique= length(stim_unique);

v_trace=cell([num_stim_unique 1]);
n_grid=length(current_template);
for i_stim = 1:num_stim_unique
    temp_trace = zeros([n_grid 1]);
    stims_temp=current_template*stim_unique(i_stim);
    temp_trace(1)=0;
    for i_t = 2:n_grid
        temp1=reshape(stims_temp(1:(i_t-1)), [i_t-1,1]);
        temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*cell_params.g), [1,i_t-1]);
        temp_trace(i_t) = temp2*temp1;
    end
    v_trace{i_stim}=temp_trace;
    fprintf('%d voltage grid done;\n', i_stim);
end

delay_prob = zeros(delay_params.n_grid+1,1);
if delay_params.type == 1 % normal
    delay_prob = normpdf((0:delay_params.n_grid),delay_params.mean,...
        delay_params.std);
elseif delay_params.type == 2 % gamma
    shape=(delay_params.mean^2)/(delay_params.std^2);
    %scale
    scale = delay_params.mean/shape;
    delay_prob = gampdf((0:delay_params.n_grid),shape,scale);
end

% we approximate the probability with densities
delay_prob = delay_prob/sum(delay_prob);
min_delay = 0;max_delay = delay_params.n_grid;
% Evaluate the probability given the gain_grid
prob_trace=zeros([num_stim_unique 1]);
for i_stim = 1:num_stim_unique
    v_trace_one=v_trace{i_stim};
    gain_one=cell_params.gain_template;
    v_th_known_one=cell_params.v_th_known;
    [prob_first_spike_delayed] = voltage_to_prob(gain_one,  v_trace_one,...
        v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
    prob_trace(i_stim) =sum(prob_first_spike_delayed);
end

end