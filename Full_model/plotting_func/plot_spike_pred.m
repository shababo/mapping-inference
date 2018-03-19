function plot_spike_pred(result_spikes,spike_curves,gain, varargin)
%plot_spike_pred(result_spikes,spike_curves,gain,figure_handle)
if ~isempty(varargin) && ~isempty(varargin{1})
    colors = varargin{1};
else
    colors=colorcube(length(gain)+1);
end
power_ref=0;

for i = 1:length(gain)
%     i=good_neurons(ii);
%     i=8;
pred_current=(result_spikes(i).power{1}-power_ref)*gain(i);
pred_current=(unique(result_spikes(i).power{1})-power_ref)*gain(i);

stim_index=zeros(length(pred_current),1);
for i_stim = 1:length(pred_current)
    [~, stim_index(i_stim)]=min(abs(pred_current(i_stim) - spike_curves.current));
    neurons(i).pred_mean(i_stim)=spike_curves.mean(stim_index(i_stim));
    neurons(i).pred_sd(i_stim)=spike_curves.sd(stim_index(i_stim));
end

act_time=[result_spikes(i).spike_times{1}];
act_time=[result_spikes(i).spike_time_means];
nanindex=find(~isnan(act_time));
scatter(act_time(nanindex)'/20,neurons(i).pred_mean(nanindex)/20,'Marker','o','SizeData',25,...
    'MarkerFaceColor',colors(i,:), 'MarkerEdgeColor',colors(i,:),'DisplayName',num2str(i));%, 'MarkerFaceAlpha',0.8)
hold on;

% x=[0 300];y=[0 300];
% line(x, y)
% hold on;
% 
% xlabel('Recorded spike time');
% ylabel('Predicted spike time (mean)');
% title(['Prediction of spike time given power'])
% 
% 
% fig.Units = 'inches';
% fig.Position = [0 0 8 4];
% 
% if ~isempty(save_path)
%     saveas(figure_index,fullfile(save_path, ['Fits' num2str(figure_index) '.png']));
%     close(figure_index)
% end


end
% legend('show')

%
% xlim([0 4]);
% ylim([0 4]);
x=[0 15];y=[0 15];
line(x, y)
hold on;

xlabel('Estimated Mean Spike Time (msec)');
ylabel('Predicted  Mean Spike Time (msec)');
title(['Prediction of Mean Spike Time Given Power and Gain (N = ' num2str(length(gain)) ', on cell target)'])


fig.Units = 'inches';
fig.Position = [0 0 8 4];

