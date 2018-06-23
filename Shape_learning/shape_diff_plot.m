function [figure_handle]=shape_diff_plot(shifts, neurons,fixed_grid, n_curves,...
    figure_handle, tau,  sigma,...
    mean_params,var_params,ax)


eig_epsilon=1e-5;
n_cell=length(neurons);
colors=lines(n_cell);
corrected_grid=cell([n_cell 1]);
for i_cell = 1:n_cell
    corrected_grid{i_cell}=[neurons(i_cell).stim_grid - shifts(i_cell)];
    X=[corrected_grid{i_cell}]';
    Y=[neurons(i_cell).scaled_current]';
    
    Xstar=fixed_grid';
    
    
    % Correlation between new points and the old points:
    nsq=sum(X.^2,2);
    nsqstar=sum(Xstar.^2,2);
    K_pred = nsqstar*ones(1,length(nsq)) +  ones(length(nsqstar),1)*nsq';
    K_pred=bsxfun(@minus,K_pred,(2*Xstar)*X.');
    
    sigma_varstar= sqrt(max(eig_epsilon,quick_match(Xstar,var_params)));
    sigma_var= sqrt(max(eig_epsilon,quick_match(X,var_params)));
    
    sigma_mat_star_left= sigma_varstar*ones(1,length(nsq));
    sigma_mat_star_right= ones(length(nsqstar),1)*sigma_var';
    
    
    K=bsxfun(@plus,nsq,nsq');
    K=bsxfun(@minus,K,(2*X)*X.');
    sigma_mat = sigma_var*ones(1,length(X));
    K=sigma_mat.*exp(-K/tau).*sigma_mat';
    K_exp=K+ diag(ones(length(X),1))*sigma(i_cell);

    K_pred= sigma_mat_star_left.*exp(-K_pred/tau).*sigma_mat_star_right;
    
    post_mean= K*inv(K_exp)*(Y-quick_match(X,mean_params))+quick_match(X,mean_params);
    pred_mean= K_pred*inv(K_exp)*(Y-quick_match(X,mean_params))+quick_match(Xstar,mean_params);
    
    % scatter(X,post_mean)
    [uniq_grid, ix]=sort(X);
    
    plot(fixed_grid,pred_mean,'Color',colors(i_cell,:))
    hold on;
    scatter(neurons(i_cell).stim_grid- shifts(i_cell),neurons(i_cell).scaled_current,...
        'MarkerFaceColor',colors(i_cell,:),'MarkerEdgeColor',colors(i_cell,:));
    hold on;
    
end
      xlim([min(fixed_grid) max(fixed_grid)])
      ylim([0 1.5])

xlabel('Stimulation location (adjusted)')
ylabel('Scaled current')
title([ax ' ' 'GP: Tau ' num2str(round(sqrt(tau),1)) '; Mean Noise sigma: ' ...
    num2str(round(mean(sigma),3))])



