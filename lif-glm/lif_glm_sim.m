function [V_vect, spikes] = lif_glm_sim(stim,params,funcs)






last_spike = 1;
dt = 1;
t_end = size(stim,2);
t_vect=0:dt:t_end;
V_spike = 70;  
spikes = zeros(size(t_vect));
V_vect=zeros(1,length(t_vect));
V_plot_vect=zeros(1,length(t_vect));
lambda = zeros(1,length(t_vect));
i=1; %index denoting which element of V is being assigned

V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V
NumSpikes=0; %holds number of spikes that have occurred

tao=exprnd(1);
%     lambda(i)=log(exp(V_vect(i)-V_th_sim)+1);
lambda(i) = exp(V_vect(i)-params.V_th);

for t=dt:dt:t_end %loop through values of t in steps of df ms        
    %V_inf = E_L + I_e_vect(i)*R_m;

    V_vect(i+1) = V_vect(i) + ...
        -params.g*V_vect(i) + stim(i)*params.gain; %Euler's method
%         lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
    lambda(i+1) = funcs.invlink(V_vect(i)-params.V_th);

    %if statement below says what to do if voltage crosses threshold
    if sum(lambda(last_spike+1:i+1))>tao %cell spiked
%         if rand < lambda(i+1)
        V_vect(i+1)=params.V_reset; %set voltage back to V_reset_sim
        V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here

        NumSpikes=NumSpikes+1; %add 1 to the total spike count
        spikes(i)=1;
        if last_spike == i+1
            disp('two spikes in a row!!')
            return
        end
        last_spike = i+1;
        tao=exprnd(1);
%             lambda(1:i)=0;
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
    else %voltage didn't cross threshold so cell does not spike
        V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th_sim)+1);
    end
    i=i+1; %add 1 to index,corresponding to moving forward 1 time step
end    

