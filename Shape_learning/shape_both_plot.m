function [figure_handle]=shape_both_plot(neurons,fixed_grid, n_curves,...
    figure_handle, tau, gains,...
    mean_params,var_params,ax)

eig_epsilon=1e-5;
subplot(1,2,1)  
n_cell=length(neurons);
colors=lines(n_cell);
corrected_grid=cell([n_cell 1]);

% check if all gains are equal:
if var(gains)==0
    fit_gain = false;
else
    fit_gain = true;
end
max_Y=0;
for i_cell = 1:n_cell
    corrected_grid{i_cell}=neurons(i_cell).adjusted_grid;
            
    X=[corrected_grid{i_cell}]';
    Y=[neurons(i_cell).scaled_current]';
    
    Xstar=fixed_grid';
    
    white_sigma = neurons(i_cell).noise_sigma;
    
     % Correlation between new points and the old points:
    nsq=sum(X.^2,2);
    nsqstar=sum(Xstar.^2,2);
    K_pred = nsqstar*ones(1,length(nsq)) +  ones(length(nsqstar),1)*nsq';
    K_pred=bsxfun(@minus,K_pred,(2*Xstar)*X.');
    
    sigma_varstar= sqrt(max(eig_epsilon,quick_match(Xstar,var_params)));
    sigma_var= sqrt(max(eig_epsilon,quick_match(X,var_params)));
     Ymean = quick_match(X,mean_params);
    if fit_gain
        
        Y=neurons(i_cell).raw_current';
         sigma_varstar =gains(i_cell)* sigma_varstar;
         sigma_var =gains(i_cell)* sigma_var;
         white_sigma = gains(i_cell)*white_sigma;
         Ymean=Ymean*gains(i_cell);
    end
    
    sigma_mat_star_left= sigma_varstar*ones(1,length(nsq));
    sigma_mat_star_right= ones(length(nsqstar),1)*sigma_var';
    
    
    K=bsxfun(@plus,nsq,nsq');
    K=bsxfun(@minus,K,(2*X)*X.');
    sigma_mat = sigma_var*ones(1,length(X));
    K=sigma_mat.*exp(-K/tau).*sigma_mat';
    K_exp=K+diag(neurons(i_cell).noise_sigma.^2);

    K_pred= sigma_mat_star_left.*exp(-K_pred/tau).*sigma_mat_star_right;
    
    post_mean= K*inv(K_exp)*(Y-Ymean)+quick_match(X,mean_params);
    pred_mean= K_pred*inv(K_exp)*(Y-Ymean)+quick_match(Xstar,mean_params);
    
    % scatter(X,post_mean)
    [uniq_grid, ix]=sort(X);
    
    plot(fixed_grid,pred_mean,'Color',colors(i_cell,:))
    hold on;
    scatter(neurons(i_cell).adjusted_grid,Y,...
        'MarkerFaceColor',colors(i_cell,:),'MarkerEdgeColor',colors(i_cell,:));
    hold on;
    max_Y=max([max_Y; Y]);
    
end
% xlim([min(fixed_grid) max(fixed_grid)])
xlim([-150 150])
ylim([0 max_Y])

xlabel('Stimulation location (adjusted)')
ylabel('Scaled current')
title([ax ' ' 'GP: Tau ' num2str(round(sqrt(tau),1)) '; Mean Shift: ' num2str(round(mean([neurons(:).initial_shift] ),1) ) ])

subplot(1,2,2) 
colors=lines(n_curves);
for i_cell = 1:n_curves

    X=fixed_grid';
    nsq=sum(X.^2,2);
    
    K=bsxfun(@plus,nsq,nsq');
    K=bsxfun(@minus,K,(2*X)*X.');
    sigma_var= sqrt(max(eig_epsilon,quick_match(X,var_params)));
    sigma_mat = sigma_var*ones(1,length(X));
    K=sigma_mat.*exp(-K/tau).*sigma_mat';

    Y =mvnrnd(quick_match(X,mean_params),K);
    plot(fixed_grid,Y,'Color',colors(i_cell,:))
    hold on;
%     if i_cell == 1
%         plot(fixed_grid,Y,'Color',colors(i_cell,:),'LineWidth',2)
%         hold on;
%     Y_obs =Y+normrnd(0, mean(sigma),[length(Y) 1]);
%         scatter(fixed_grid, Y_obs,...
%         'MarkerFaceColor',colors(i_cell,:),'MarkerEdgeColor',colors(i_cell,:));
%     hold on;
%     
%     end
end    


Y =quick_match(X,mean_params);
    plot(fixed_grid,Y,'Color','k','LineWidth',3)
    hold on;
%     xlim([min(fixed_grid) max(fixed_grid)])
      ylim([0 1.5])
% xlim([-30 30])
xlim([-150 150])
xlabel('Stimulation location')
ylabel('Current')
title(['Simulated (' ax ')'])
