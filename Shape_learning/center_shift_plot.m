function [figure_handle]=center_shift_plot(true_shift, est_shift, figure_handle,method)


mse= (mean((true_shift-est_shift).^2));
scatter(true_shift,est_shift,'MarkerFaceColor','r')
hold on;
line([-100 100], [-100 100])
xlim([-40 40])
ylim([-40 40])
xlabel('True shifts')
ylabel('Est. shifts')
title([method ' ' num2str(round(mse,0))])