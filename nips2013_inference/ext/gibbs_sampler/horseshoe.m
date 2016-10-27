function [Wt,W2t,identifier]=horseshoe(M,r,W,nG,sample_sigma,sample_tau,sigma,tau)
% HORSESHOE a Gibbs sampler of horseshoe model as described in Scott's 2009 thesis. 
% Implements Ari's idea of the alternating envelopes in the rejection sampling.

%% Set flags
plot_estimate=0; % want to plot the estimate?
plot_sample_cloud=0; % want to plot a cloud of projected W samples?
plot_sample_trajectories=0; % want to plot Gibbs trajectories of the first few weight components?

if nargin<3 error('horseshoe:too_few_arguments','not enough input arguments'); end
if nargin<4 nG=1000; end % number of Gibbs sweeps after burn-in
if nargin<5 sample_sigma=0; end
if nargin<6 sample_tau=0; end
if nargin<7 sigma=1; end
if nargin<8 tau=1; end

%% Initialize Gibbs sampler
p=length(W);
lambda=abs(tau*tan(pi*(rand(p,1)-0.5))); % sample from half-Cauchy prior
% lambda=W+randn(p,1); % start with this value for lambda as a cheat to
% avoid waiting a long time for the sampler to find a region of high
% posterior density with high acceptance rates for rejection sampling of
% lambda. i get similar results as from starting with a draw from the prior
% over lambda, except that that way takes longer at first.

%% Estimate \W_i from y_i by Gibbs sampling
Wt=zeros(p,1); % W tally
W2t=zeros(p,1); % W.^2 tally
nB=nG/10; % # of burn-in sweeps
if (plot_sample_cloud || plot_sample_trajectories) Ws=zeros(p,nG); end % array of samples
for g=-nB+1:nG
    disp(['Gibbs sweep ' num2str(g) ' of ' num2str(nG)])

    A=M/sigma^2+1/tau^2*diag(lambda.^-2);
    A_inv_r=(A\r)/sigma^2;
    Wg=mvnrnd(A_inv_r,inv(A))';
    if g>0 % after burn-in, start recording
        Wt=Wt+A_inv_r;
        W2t=W2t+Wg.^2;
        if (plot_sample_cloud || plot_sample_trajectories) Ws(:,g)=Wg; end
    end
    
    for j=1:p
        b=Wg(j)^2/(2*tau^2);
        L=full_conditional_sampler_lambda(b);
        lambda(j)=sqrt(L);
    end
    
    if sample_tau
        a=(p-1)/2; % a parameter of gamma-like distribution
        b=sum(Wg.^2./lambda.^2)/(2*sigma^2);
        L=full_conditional_sampler_tau(a,b);
        tau=sqrt(L);
    end
    
    if sample_sigma
        a=(p-1)/2;
        b=1;
        L=full_conditional_sampler_sigma(a,b);
        sigma=tau/sqrt(L);
    end
end
Wt=Wt/nG;
W2t=W2t/nG-Wt.^2;

%% Make labels to identify printed figures
[~,timestamp]=system(['echo -n `date +day_%D_time_%T | sed ''s/[\:\/]/\_/g''`']);
identifier=['nG' num2str(nG)];

%% Graphical output
if plot_estimate feval(plotting_horseshoe('plot_estimate'),W,Wt,W2t); end
if plot_sample_cloud feval(plotting_horseshoe('plot_sample_cloud'),W,Wt,W2t,Ws); end
if plot_sample_trajectories feval(plotting_horseshoe('plot_sample_trajectories'),Ws); end

end