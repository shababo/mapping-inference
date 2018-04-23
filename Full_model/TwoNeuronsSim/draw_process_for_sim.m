function [mpp]=draw_process_for_sim(mpp,neurons,spike_curves)
Ntrials = length(mpp);
bg_rate=1e-5;
Tmax=300;
for i= 1:Ntrials
   spike_times = [];
   event_times = [];
   assignments=[];
   eff_size=[];
   for i_neuron = 1:length(neurons)
      actual_stim= neurons(i_neuron).gain*mpp(i).true_stimulation(i_neuron);
      [~, Ia]=min(abs(actual_stim - spike_curves.current));
      spike_param = struct;
      spike_param.mean=spike_curves.mean(Ia);
      spike_param.sd=spike_curves.sd(Ia);
      spike_one = normrnd(spike_param.mean,spike_param.sd);
      delay_one = normrnd(neurons(i_neuron).delay_mean,sqrt(neurons(i_neuron).delay_var));
      % truncate the event if it is larger than Tmax
      if ((spike_one+delay_one)<Tmax) & (unifrnd(0,1)<neurons(i_neuron).PR)
          spike_times = [spike_times spike_one];
          assignments=[assignments i_neuron];
          event_times = [event_times spike_one+delay_one];
          eff_size=[eff_size mpp(i).stimulation(i_neuron)];
      end
   end
   % injecting background events:
   bg_prob=bg_rate*Tmax;
   bg_yes = binornd(1,bg_prob);
   if bg_yes
       spike_one=unifrnd(0,Tmax);
       spike_times = [spike_times spike_one];
       assignments=[assignments 0];
       event_times = [event_times spike_one];
       eff_size=[eff_size 0];
   end
   
   mpp(i).event_times = event_times;
   mpp(i).spike_times = spike_times;
   mpp(i).assignments = assignments;    
   mpp(i).eff_size = eff_size;    
end