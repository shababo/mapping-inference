function [figure_handle]=shape_sim_plot(fixed_grid, n_curves, figure_handle,method,tau, sigma_gp, sigma,...
    mean_func,ax)
 

colors=lines(n_curves);
for i_cell = 1:n_curves

    X=fixed_grid';
    
    nsq=sum(X.^2,2);
    K=bsxfun(@plus,nsq,nsq');
    K=bsxfun(@minus,K,(2*X)*X.');
    K=sigma_gp*exp(-K/tau);
    
    Y =mvnrnd(mean_func(X),K);
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


Y =mean_func(X);
    plot(fixed_grid,Y,'Color','k','LineWidth',3)
    hold on;
    
      
% xlim([-30 30])

xlabel('Stimulation location')
ylabel('Current')
title(['Simulated (' ax ') ' 'GP: Tau ' num2str(round(tau,1)) ', ' 'Sigma ' num2str(round(sigma_gp,3)) '; Mean Noise sigma: ' ...
    num2str(round(mean(sigma),3))])

