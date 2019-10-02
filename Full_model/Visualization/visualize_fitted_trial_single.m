function [covered_flag] = visualize_fitted_trial_single(this_trial, plot_params)
alpha=plot_params.alpha;
if plot_params.prediction
    chosen_field = 'predicted';
else
      chosen_field = 'fitted';
end
covered_flag=1;
realtimepoints=this_trial.(chosen_field).timepoints/this_trial.(chosen_field).time_factor;
n_neuron= size(this_trial.(chosen_field).intensity.event,1); % This will always include the background rate 
neuron_colors = lines(n_neuron);

if plot_params.spike 
    fits = this_trial.(chosen_field).intensity.spike;
else
    fits =this_trial.(chosen_field).intensity.event;
end

if plot_params.stack 
    % Stack the intensities 
    fits = sum(fits);
end
if isfield(plot_params, 'LineWidth')
    lw=plot_params.LineWidth;
else
    lw=5;
end
if isfield(plot_params,'shift')
    shift=plot_params.shift;
else
    shift=0;
end 
if    isfield(plot_params,'vertical_gap') 
    barheight=plot_params.gap;
else
    barheight = plot_params.max_intensity;
end

if  isfield(plot_params,'max_intensity') 
    max_int=plot_params.max_intensity;
else
    max_int = max(max(fits));
end


if isfield(plot_params, 'quantiles')
    quantile_prob=plot_params.quantiles;
else
    quantile_prob = [0.25 0.75];
end

for i = 1:size(fits,1)
    switch plot_params.fit_type
        case 'full_intensity'
            total_int = cumsum(fits(i,:));
            total_int=total_int/max(total_int);
            plot(fits(i,:)/max_int*plot_params.gap+shift*ones(1,size(fits,2)),realtimepoints,'Color',plot_params.colors(plot_params.this_loc,:),'LineWidth',lw)
            qt_lower=max(realtimepoints(find(total_int<quantile_prob(1)) ));
              if isempty(qt_lower)
                  qt_lower=min(realtimepoints);
              end
              qt_upper=max(realtimepoints(find(total_int<quantile_prob(2)) ));
              if isempty(qt_upper)
                  qt_upper=max(realtimepoints);
              end
              quantiles_tmp =[qt_lower qt_upper];
            scatter([shift  shift], quantiles_tmp, plot_params.markerSize,  '^','MarkerFaceColor','red',...
                'MarkerEdgeColor','red')
        
        case 'quantiles'
            total_int = cumsum(fits(i,:));
              total_int=total_int/max(total_int);
            
              qt_lower=max(realtimepoints(find(total_int<quantile_prob(1)) ));
              if isempty(qt_lower)
                  qt_lower=min(realtimepoints);
              end
              qt_upper=max(realtimepoints(find(total_int<quantile_prob(2)) ));
              if isempty(qt_upper)
                  qt_upper=max(realtimepoints);
              end
              quantiles_tmp =[qt_lower qt_upper];
              plot([shift  shift], quantiles_tmp,':','LineWidth',lw+1,'Color',plot_params.colors(plot_params.this_loc,:))
           
    end
    
    hold on;
end

if isfield(this_trial,'event_times')
    if ~isempty(this_trial.event_times)
        tmp=this_trial.event_times/this_trial.(chosen_field).time_factor;
        scatter(shift, tmp,'MarkerFaceColor',plot_params.colors(plot_params.this_loc,:),...
              'MarkerEdgeColor',plot_params.colors(plot_params.this_loc,:))
        if strcmp(plot_params.fit_type,'quantiles')
            covered_flag = (quantiles_tmp(1) < tmp ) & (quantiles_tmp(2) > tmp );
        end
        
        
        %plot([tmp;tmp], [vertical_shift;vertical_shift+barheight],'-','LineWidth',lw+1,'Color',plot_params.colors(plot_params.this_loc,:))
        hold on;
    end

end


if (~covered_flag) & strcmp(plot_params.fit_type,'quantiles')
    text(shift,tmp+0.2,['Loc ' num2str( plot_params.this_loc) ' Pow ' num2str(this_trial.power_levels)]);
end

if strcmp(plot_params.fit_type,'quantiles')
plot([shift  shift], [0 sum(sum(fits))],'-','LineWidth',lw+5,'Color',plot_params.colors(plot_params.this_loc,:))
end
% if ~isfield(plot_params,'shift')
%     xlabel('Time (ms)');
%     if plot_params.prediction
%     ylabel('Predicted intensity')
%     else
%         ylabel('Fitted intensity')
%     end
% end

   
% Visualize multiple trials at the same locations, ordered by power levels 