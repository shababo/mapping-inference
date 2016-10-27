%% Input parameters
get_new_data=0; % want to generate new M and r?
save_this_data=1; % want to save this M and r for the future? (will overwrite previous M and r)
run_horseshoe=0; % want to run the horseshoe sampler?
    hs_sample_sigma=0; % want to sample sigma in the horseshoe sampler?
    hs_sigma=1; % if not, what value is sigma fixed to? (irrelevant if not sampling sigma)
    hs_sample_tau=0; % want to sample tau in the horseshoe sampler?
    hs_tau=1; % if not, what value is tau fixed to? (irrelevant if not sampling tau)
run_spike_and_slab=1; % want to run the spike-and-slab sampler?
    ss_alpha_a=1; % alpha parameter of beta prior on 'a'
    ss_beta_a=1; % beta parameter of beta prior on 'a'
    ss_alpha_tau=1; % alpha parameter of inverse gamma prior on tau
    ss_beta_tau=1; % beta parameter of inverse gamma prior on tau
run_spike_and_slab_exponential=0;
    sse_lambda=30;
plot_horseshoe=0; % want to plot the output of the horseshoe sampler?
plot_spike_and_slab=1; % want to plot the output of the spike-and-slab sampler?
print_output=1; % want to print the figure to pdf?
nG=100; % # of Gibbs sweeps after burn-in period

%% Main program
if get_new_data sampler_run;
else load('D'); end
if save_this_data save('D','M','r','true_wts'); end
W=true_wts; % W is the true weights vector
if run_horseshoe [Wt,W2t,id_HS]=horseshoe(M,r,W,nG,hs_sample_sigma,hs_sample_tau,hs_sigma,hs_tau); end
if run_spike_and_slab [Qs,id_SS]=spike_and_slab(M,r,W,nG,ss_alpha_a,ss_beta_a,ss_alpha_tau,ss_beta_tau); end
if run_spike_and_slab_exponential [Qs,id_SS]=spike_and_slab_exponential(M,r,W,nG,ss_alpha_a,ss_beta_a,sse_lambda); end

%% Graphical output
params=struct;
if plot_horseshoe 
    params.ssig=hs_sample_sigma;
    params.sig=hs_sigma;
    params.stau=hs_sample_tau;
    params.tau=hs_tau;
    params.Wt=Wt;
    params.W2t=W2t;
end
if plot_spike_and_slab
    params.Qs=Qs;
    params.aa=ss_alpha_a;
    params.ba=ss_beta_a;
    params.at=ss_alpha_tau;
    params.bt=ss_beta_tau;
end
if (plot_horseshoe || plot_spike_and_slab)
    feval(plotting_horseshoe('plot_horseshoe_and_spike_and_slab'), ...
        plot_horseshoe,plot_spike_and_slab,print_output,W,M,r,params);
end