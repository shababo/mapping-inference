function plot_event_pred(delay_params,spike_curves,varargin)
%%
%     ax1=axes('Position',[0.08 0.08 0.85 0.85],'Visible','off');
%     axes(ax1)
    avai_index=[];
    for i = 1:length(delay_params.mpp)
        if ~isnan(delay_params.mpp(i).event_times)
        avai_index = [avai_index i];
        end
    end
    
    event_times=[delay_params.mpp(avai_index).event_times];
    stim_size=[delay_params.mpp(avai_index).stimulation];
    
    scatter(stim_size,event_times/20,'o','SizeData',20,...
        'MarkerFaceColor','b',...
        'MarkerEdgeColor','b',...
        'MarkerFaceAlpha',1)
    xlim([0 100])
    ylim([0 30])
    hold on;
    
    %%
    % Mean of spike time:
    plot(spike_curves.current/delay_params.gain,spike_curves.mean/20,'Color','r','LineWidth',2)
    hold on;
    % Mean of event time:
    predicted_mean=spike_curves.mean+delay_params.mean;
    plot(spike_curves.current/delay_params.gain,predicted_mean/20,'Color','b','LineWidth',2)
    hold on;
    %%
    % The error bar (spike var)
    sd_spike=sqrt(spike_curves.sd.^2);
    spike_upp=spike_curves.mean+sd_spike;
    plot(spike_curves.current/delay_params.gain,spike_upp/20,'Color','r','LineWidth',1,...
        'LineStyle',':')
    hold on;
    spike_low=spike_curves.mean-sd_spike;
    plot(spike_curves.current/delay_params.gain,spike_low/20,'Color','r','LineWidth',1,...
        'LineStyle',':')
    hold on;
    
    sd_event=sqrt(spike_curves.sd.^2+delay_params.var);
    event_upp=predicted_mean+sd_event;
    plot(spike_curves.current/delay_params.gain,event_upp/20,'Color','b','LineWidth',1,...
        'LineStyle',':')
    hold on;
    event_low=predicted_mean-1*sd_event;
    index_low=find(event_low<1);
    event_low(index_low)=1;
    plot(spike_curves.current/delay_params.gain,event_low/20,'Color','b','LineWidth',1,...
        'LineStyle',':')
    hold on;
    %%
    xlabel('Power');
    ylabel('PSC event time');
    
    