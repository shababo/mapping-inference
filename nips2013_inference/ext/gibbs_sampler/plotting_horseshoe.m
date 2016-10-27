function h_pf=plotting_horseshoe(f_str)
% PLOTTING_HORSESHOE plotting functions for horseshoe

if strcmp(f_str,'plot_estimate')
    h_pf=@plot_estimation_of_first_few;
elseif strcmp(f_str,'plot_sample_cloud')
    h_pf=@plot_sample_cloud;
elseif strcmp(f_str,'plot_sample_trajectories')
    h_pf=@plot_sample_trajectories;
elseif strcmp(f_str,'plot_horseshoe_and_spike_and_slab')
    h_pf=@plot_horseshoe_and_spike_and_slab;
else
    save('error')
    error(['plotting_horseshoe: requested function ' f_str ' does not exist'])
end

end


function plot_estimate(theta,thetat,theta2t,y)
% PLOT_ESTIMATE compares actual, observed, and estimated values theta.

fig_params=initalize_fig_params();
h_fig=figure('name','Estimate','numbertitle','off'); clf; hold on
prepare_plot(h_fig,13,5);
plot(find(theta~=0),theta(theta~=0),'linestyle','none','marker','o','markersize',15,'color','r','displayname','$W$','linewidth',1)
if nargin>3 plot(y,'linestyle','none','marker','x','markersize',20,'color','k','displayname','$y$','linewidth',1); end
errorbar(thetat,2*sqrt(theta2t),'linestyle','none','marker','+','markersize',20,'color','k','displayname','$\hat{W}$','linewidth',1)
hline(0,'k--')
set(gca,'fontsize',16)
xlabel('Element $i$','fontsize',fig_params(1).label_font_size)
ylabel('$W_i$','fontsize',fig_params(1).label_font_size)
h_leg=legend('show');
set(h_leg,'fontsize',fig_params(1).legend_font_size);
set(h_leg,'interpreter','latex')
h_title=title('Horseshoe estimate: Rao-Blackwell Gibbs sampler mean','fontsize',16);
set(h_title,'fontsize',fig_params(1).title_font_size)

end


function plot_sample_cloud(theta,thetat,theta2t,thetas,y)
% PLOT_SAMPLE_CLOUD scatter plots the samples throughout the Gibbs sampler.

fig_params=initalize_fig_params();
h_fig=figure('name','Sample cloud','numbertitle','off'); clf; hold on
prepare_plot(h_fig,13,5);
set(0,'defaultaxesfontname','cmu serif')
nG=size(thetas,2);
scatter_step=floor(nG/1000);
for j=1:scatter_step:nG
    scatter(thetas(1,j),thetas(2,j),10,[0 (j-1)/nG 0],'filled')
%     drawnow
end
h_scat_1=scatter(theta(1),theta(2),75,[1 0 0],'filled','displayname','$\theta$');
h_scat_2=scatter(thetat(1),thetat(2),75,[0 0 1],'filled','displayname','$\hat{\theta}$');
if nargin>4 h_scat_3=scatter(y(1),y(2),75,[0 0 0],'filled','displayname','$y$'); end
h_leg=legend([h_scat_1 h_scat_2]);
set(h_leg,'interpreter','latex')
h_err=herrorbar(thetat(1),thetat(2),2*sqrt(theta2t(1)));
set(h_err,'color','b','linewidth',1)
errorbar(thetat(1),thetat(2),2*sqrt(theta2t(2)),'linestyle','none','linewidth',1,'color','b')
hline(0,'k--')
vline(0,'k--')
% min_x=min(min(theta(1),thetat(1)),y(1));
% max_x=max(max(theta(1),thetat(1)),y(1));
% min_y=min(min(theta(2),thetat(2)),y(2));
% max_y=max(max(theta(2),thetat(2)),y(2));
% axis([min_x-(max_x-min_x)/2 max_x+(max_x-min_x)/2 min_y-(max_y-min_y)/2 max_y+(max_y-min_y)/2])
set(gca,'fontsize',16)
title('Close-up of marginal sampling distribution of $\theta_1$ and $\theta_2$','fontsize',16)

end


function plot_sample_trajectories(thetas)
% PLOT_SAMPLE_TRAJECTORIES plots the sample trajectories of the Gibbs 
% samples of first four elements of $\theta$.

fig_params=initalize_fig_params();
h_fig=figure('name','Sample trajectories','numbertitle','off'); clf; hold on
prepare_plot(h_fig,13,5);
set(0,'defaultaxesfontname','cmu serif')
plot(thetas(1,:),'r','displayname','$\hat{\theta}_1$')
plot(thetas(2,:),'k','displayname','$\hat{\theta}_2$')
plot(thetas(3,:),'b','displayname','$\hat{\theta}_3$')
plot(thetas(4,:),'g','displayname','$\hat{\theta}_4$')
h_title=title('The sampling paths of the components of $\theta$');
set(h_title,'fontsize',fig_params(1).title_font_size)
set(gca,'fontsize',16)
h_leg=legend('show');
set(h_leg,'interpreter','latex')

end


function plot_horseshoe_and_spike_and_slab(phs,pss,po,W,M,r,pars)
% PLOT_HORSESHOE_AND_SPIKE_AND_SLAB plots estimates of weights for
% horseshoe or spike-and-slab or both against ML estimate and true weights.

h_fig=figure('name','HS SS ML estimates','numbertitle','off'); clf
prepare_plot(h_fig,8.5,5.5*(phs+pss));
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','cmu serif')

if phs && pss subplot(2,1,1); end
if phs
    hold on
    plot(find(W~=0),W(W~=0),'r.','linestyle','none','markersize',45,'displayname','$W$')
    plot(M\r,'g.','displayname','$\hat{W}_{ML}$','markersize',45)
    errorbar(pars.Wt,2*sqrt(pars.W2t),'k.','linestyle','none','markersize',30,'displayname','$\hat{W}_{HS}$')
    set(gca,'fontsize',16)
    if pars.ssig sigma_str='$\sigma$ sampled';
    else sigma_str=['$\sigma=' num2str(pars.sig) '$']; end
    if pars.stau tau_str='$\tau$ sampled';
    else tau_str=['$\tau=' num2str(pars.tau) '$']; end
    title(['Horseshoe, ' sigma_str ', ' tau_str],'fontsize',16);
    if ~pss xlabel('Element $i$','fontsize',16); end
    ylabel('$W_i$','fontsize',16)
    h_leg=legend('show');
    set(h_leg,'interpreter','latex','fontsize',16)
    for k=1:3 hline(k,'k:'); end
    hline(0,'k--');
    if any(pars.Wt==0) plot(find(pars.Wt==0),0,'w.','markersize',15); end
    axis([0 36 -1 4])
end
if (phs && pss) subplot(2,1,2); end
if pss
    hold on
    plot(find(W~=0),W(W~=0),'r.','linestyle','none','markersize',45)
    plot(M\r,'g.','displayname','$\hat{W}_{ML}$','markersize',45)
    errorbar(1:size(pars.Qs,1),pars.Qs(:,2),pars.Qs(:,2)-pars.Qs(:,1), ...
        pars.Qs(:,3)-pars.Qs(:,2),'k.','linestyle','none','markersize',30)
    set(gca,'fontsize',16)
    a_str=['$\alpha_a=' num2str(pars.aa) '$, $\beta_a=' num2str(pars.ba) '$'];
    tau_str=['$\alpha_\tau=' num2str(pars.at) '$, $\beta_\tau=' num2str(pars.bt) '$'];
    title(['Spike-and-slab, ' a_str ', ' tau_str],'fontsize',16)
    xlabel('Element $i$','fontsize',16)
    ylabel('$W_i$','fontsize',16)
    for k=1:3 hline(k,'k:'); end
    hline(0,'k--');
    if any(pars.Qs(:,2)==0) plot(find(pars.Qs(:,2)==0),0,'w.','markersize',15); end
    axis([0 36 -1 4])
end
directory='~/research/dendritic_synaptic_connectivity/matlab/figures/compare_horseshoe_spike_and_slab';
identifier='';
if phs
    if pars.ssig identifier=[identifier 'HS_sX_'];
    else identifier=[identifier 'HS_s' num2str(pars.sig) '_'];
    end
    if pars.stau identifier=[identifier 'tX_'];
    else identifier=[identifier 't' num2str(pars.tau)];
    end
end
if (phs && pss) identifier=[identifier '_']; end
if pss identifier=[identifier 'SS_aa' num2str(pars.aa) '_ba' num2str(pars.ba) ...
        '_at' num2str(pars.at) '_bt' num2str(pars.bt)]; end
[~,timestamp]=system('echo -n `date +D%D_T%T | sed ''s/[\:\/]//g''`');
% if print_output laprint(h_fig,[directory '/' identifier '_' timestamp]); end
% if po print('-dpdf',[directory '/' identifier '_' timestamp]); end

end


function prepare_plot(h_fig,width,height)
% PREPARE_PLOT prints the current plot to a pdf file

if nargin<2 width=10; end
if nargin<3 height=6; end
orient portrait
set(h_fig,'units','inches')
set(h_fig,'outerposition',[0 0 width height])
set(h_fig,'paperpositionmode','auto');
p=get(h_fig,'position');
set(h_fig,'PaperSize',[p(3) p(4)]);

end


function fig_params=initalize_fig_params()
% FIG_PARAMS sets default values for figure parameters such as font sizes

fig_params=struct('axis_font_size',12,'label_font_size',18, ...
    'title_font_size',20,'colors',{'k','r','g','b'}, ...
    'legend_font_size',16);

end


function latex_ticks(h_ax)
% LATEX_TICKS make tick labels LaTeX

xs_old=get(h_ax,'xtick');
xs_new={};
h_lab_x=get(h_ax,'xlabel');
pos_x=get(h_lab_x,'position');
for (j=1:length(xs_old)) xs_new{j}=['$' num2str(xs_old(j)) '$']; end
ys_old=get(h_ax,'ytick');
ys_new={};
h_lab_y=get(h_ax,'ylabel');
pos_y=get(h_lab_y,'position');
for (j=1:length(ys_old)) ys_new{j}=['$' num2str(ys_old(j)) '$']; end
format_ticks(h_ax,xs_new,ys_new);
set(h_lab_x,'position',pos_x)
set(h_lab_y,'position',pos_y)

end

