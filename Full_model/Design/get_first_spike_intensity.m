function [prob_trace,v_trace] = get_first_spike_intensity(...
    linkfunc,...
    current_template,stim_grid,cell_params,delay_params)
%----------- Parameters of the template cell

g=cell_params.g;

n_grid=length(current_template);
num_stim_grid= length(stim_grid);
v_trace=cell([num_stim_grid 1]);
for i_stim = 1:num_stim_grid
    temp_trace = zeros([n_grid 1]);
    stims_temp=current_template*stim_grid(i_stim);
    temp_trace(1)=0;
    for i_t = 2:n_grid
        temp1=reshape(stims_temp(1:(i_t-1)), [i_t-1,1]);
        temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*g), [1,i_t-1]);
        temp_trace(i_t) = temp2*temp1;
    end
    v_trace{i_stim}=temp_trace;
%     fprintf('%d voltage grid done;\n', i_stim);
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
prob_trace=zeros(num_stim_grid,n_grid);

for i_stim = 1:num_stim_grid
    v_trace_one=v_trace{i_stim};
    optical_gain=cell_params.optical_gain;
    v_thresh=cell_params.v_thresh;
    [prob_first_spike_delayed] = voltage_to_prob(optical_gain,  v_trace_one,...
        v_thresh,linkfunc,delay_params,min_delay,max_delay,delay_prob);
    prob_trace(i_stim,:) =prob_first_spike_delayed;
end

end



