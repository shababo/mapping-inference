function [mpp]=design_trials_for_sim(neurons,loc_target,powers,weights)

Ntrials = length(powers)*length(loc_target);
mpp(Ntrials)=struct;


for i_pow=1:length(powers)
    sel_i = randsample(1:length(loc_target),length(loc_target),true,weights);
    for i_loc = 1:length(loc_target)
        
        n= (i_loc-1)*length(powers) + i_pow;
        mpp(n).power=powers(i_pow);
        mpp(n).location=sel_i(i_loc);
        mpp(n).location_coord=loc_target(mpp(n).location);
        mpp(n).true_stimulation = mpp(n).power*[neurons(1).z.true_received(mpp(n).location) neurons(2).z.true_received(mpp(n).location)];
        mpp(n).stimulation = mpp(n).power*[neurons(1).z.received(mpp(n).location) neurons(2).z.received(mpp(n).location)];
    end
end

