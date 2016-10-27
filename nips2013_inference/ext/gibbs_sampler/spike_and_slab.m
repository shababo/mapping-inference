%function [Qs,identifier]=spike_and_slab(M,r,W,nG,a_alpha,a_beta,tau_alpha,tau_beta)
function [Wgt Wgtv Sgt]=spike_and_slab(M,r,W,nG,a_alpha,a_beta,tau_alpha,tau_beta)

%% Set flags
watch=0; % whether or not to produce graphical output during the Gibbs sampler (makes it much slower but gives you an idea of what's going on)
quartiles=0; % whether or not to compute quartiles

%% Initialize sampler
p=length(W);
if nargin<4 nG=100; end % number of Gibbs sweeps after burn-in
if nargin<5 a_alpha=250; end
if nargin<6 a_beta=50; end
if nargin<7 tau_alpha=.9; end
if nargin<8 tau_beta=.9; end

%% Gibbs sampler for p(S,W|D)
nB=ceil(nG/10); % # of burn-in sweeps
a=betarnd(a_alpha,a_beta);
Sg=zeros(p,1);
Sgs=zeros(p,nG);
Wg=zeros(p,1);
Wgs=zeros(p,nG);
mus=zeros(p,nG);
Cs=zeros(p,nG);
for j=1:p
    Sg(j)=(rand>a);
end
Sgt=zeros(p,1);
Wgt=zeros(p,1);
Wg_sqt=zeros(p,1);
for g=-nB+1:nG
%     display(['Sweep number ' num2str(g) ' of ' num2str(nG)])
    % Sample $\tau$.
    tau_alpha_cond=tau_alpha+sum(Sg)/2;
    tau_beta_cond=tau_beta+1/2*sum(Wg.^2);
    tau=1/gamrnd(tau_alpha_cond,1/tau_beta_cond);
    tau_sq_inv=1/tau;
    % Sample $a$.
    a_alpha_cond=a_alpha+sum(~Sg);
    a_beta_cond=a_beta+sum(Sg);
    a=betarnd(a_alpha_cond,a_beta_cond);
    one_minus_a_over_a=(1-a)/a;
    % Sample $S$. First, propose pointwise flips
    for i=1:p
        % Sample from p(s_i|S_{-i},W_{-i})
        kg=sum(Sg)-Sg(i);
        ind=find(Sg==1);
        ind=ind(ind~=i);
        M_minus_i=M(ind,ind);
        r_minus_i=r(ind);
        m_i=M(:,i);
        m_i=m_i(ind);
        A=M_minus_i+tau_sq_inv*eye(kg);
        if kg>0
            factor_1=M(i,i)+tau_sq_inv-m_i'*(A\m_i);
            factor_2=m_i'*(A\r_minus_i)-r(i);
        else
            factor_1=M(i,i)+tau_sq_inv;
            factor_2=-r(i);
        end
        R=sqrt((2*pi)/factor_1)*exp(factor_2^2/factor_1);
        beta_i=one_minus_a_over_a*R;
        alpha_i=1/(1+beta_i);
        Sg(i)=(rand>alpha_i);
    end
    % Next, propose (nearest neighbor) pairwise swaps (to avoid getting
    % stuck or effectively stuck)
    for i=1:p-1
        if Sg(i)~=Sg(i+1)
            % s_i=0, s_{i+1}=1
            if sum(Sg)-Sg(i)>0
                ind_i=find(Sg==1);
                ind_i=ind_i(ind_i~=i);
                M_i=M(ind_i,ind_i);
                r_i=r(ind_i);
                A_i=M_i+tau_sq_inv*eye(sum(Sg)-Sg(i));
            else
                r_i=0;
                A_i=1;
            end
            % s_i=1, s_{i+1}=0
            if sum(Sg)-Sg(i+1)>0
                ind_i_plus_1=find(Sg==1);
                ind_i_plus_1=ind_i_plus_1(ind_i_plus_1~=(i+1));
                M_i_plus_1=M(ind_i_plus_1,ind_i_plus_1);
                r_i_plus_1=r(ind_i_plus_1);
                A_i_plus_1=M_i_plus_1+tau_sq_inv*eye(sum(Sg)-Sg(i+1));
            else
                r_i_plus_1=0;
                A_i_plus_1=1;
            end
            % Do the swap
            beta_pair=sqrt(det(A_i)/det(A_i_plus_1))*exp(1/2*(r_i_plus_1'*(A_i_plus_1\r_i_plus_1)-r_i'*(A_i\r_i)));
            alpha_pair=1/(1+beta_pair);
            swap=(rand>alpha_pair);
            Sg(i)=swap;
            Sg(i+1)=~swap;
        end
    end
    % Sample $W$.
    ind=find(Sg==1);
    M_ones=M(ind,ind);
    r_ones=r(ind);
    A=M_ones+tau_sq_inv*eye(sum(Sg));
    Wg=zeros(p,1);
    A_inv=inv(A);
    A_inv_r=A\r_ones;
    if sum(Sg)>0
        Wg(ind)=mvnrnd(A_inv_r,A_inv)';
    end
    if g>0
        Wgs(:,g)=Wg;
        Sgs(:,g)=Sg;
        if watch feval(plotting_spike_and_slab('plot_gibbs_samples'),W,Wgs,Y,g); end
        Sgt=Sgt+Sg;
        Wgt=Wgt+Wg;
        Wg_sqt=Wg_sqt+Wg.^2;
        v1=1;
        for j=1:p
            if Sgs(j,g)
                mus(j,g)=A_inv_r(v1);
                Cs(j,g)=A_inv(v1,v1);
                v1=v1+1;
            end
        end
    else 
        if watch feval(plotting_spike_and_slab('plot_gibbs_samples'),W,Wg,Y,g); end
    end
end
Sgt=Sgt./nG;
Wgt=Wgt./nG;
Wg_sqt=Wg_sqt/nG;
Wgtv=Wg_sqt-Wgt.^2;

% %% Graphical output
% identifier=['nG' num2str(nG)];
% [~,timestamp]=system('echo -n `date +day_%D_time_%T | sed ''s/[\:\/]/\_/g''`');
% % feval(plotting_spike_and_slab('plot_sample_autocorrelation_function'),Wgs)
% feval(plotting_spike_and_slab('plot_average_sparsity_vector'),Sgt);
% if quartiles
%     Qs=compute_quartiles(Sgs,mus,Cs);
% end
% feval(plotting_spike_and_slab('plot_weights_reconstruction'), ...
%     W,Wgt,Wgtv,identifier,timestamp);

end