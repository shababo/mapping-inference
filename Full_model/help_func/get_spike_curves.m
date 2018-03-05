function [spike_curves, x_current, y_spike_mean] = get_spike_curves(result_current,result_spikes,varargin)
% inputs are the two result_current and result_spikes data structures 
% call get_spike_curves(result_current,result_spikes)
% To plot the fitted curve, call
% figure_index=1;
% figure(figure_index)
% scatter(spike_curves.x_current,spike_curves.y_spike_mean)
% hold on;
% plot(spike_curves.current*1e3,spike_curves.mean)
% hold off;
% xlim([400 3000])
% xlabel('Actual current');ylabel('Mean spike time');
% 
x_current=[];
y_spike_mean=[];
y_spike_sd=[];
    
for i_cell = 1:length(result_current)
    if ~isempty(result_current(i_cell).peak_current_means)
    i_cell
    %     if i_cell ~= 11
%     if length(result_current(i_cell).these_powers) ==  length(result_spikes(i_cell).these_powers)
%         if min(result_current(i_cell).these_powers ==  result_spikes(i_cell).these_powers)==1 
            % NOTE: using these two criteria will leave only 4 cells.
            % It seems that not all the cells have matching power levels
            % between the two structures

        [~,i_c,i_s] = intersect(result_current(i_cell).these_powers,result_spikes(i_cell).these_powers);
            % make sure the power levels line with those in result_current and
            % result_spikes
            x_current =[ x_current result_current(i_cell).peak_current_means(i_c)];
            y_spike_mean=[y_spike_mean result_spikes(i_cell).spike_time_means(i_s)];
            y_spike_sd=[y_spike_sd result_spikes(i_cell).spike_time_jitter(i_s)];
%         end
%     end
    %     end
    end
end

if ~isempty(varargin) && ~isempty(varargin{1})
    specs= varargin{1};
else
    specs=struct;

%     specs.F_mean = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3);
    specs.F_mean = @(x,xdata) min(y_spike_mean) + x(1)./(xdata);
    specs.F_sd = @(x,xdata) min(y_spike_sd) + x(1)./(xdata);

    specs.current_multiplier=1e-3;
    specs.current_min=10;
    specs.current_max=3000;
    specs.current_gap=2;

    specs.sd_min=1;
    specs.sd_max=5;


end





%% Call fmincon:

ni_mean=isnan(y_spike_mean) | y_spike_mean > 160 | x_current > 3500 | isnan(x_current);
assignin('base','y_spike_mean',y_spike_mean)
assignin('base','x_current',x_current)
find(ni_mean)
% ni_mean=isnan(y_spike_mean) | x_current > 3500;
xdata=x_current(~ni_mean);ydata=y_spike_mean(~ni_mean);

x0=[120 .01 min(y_spike_mean)];
x0 = [1 1];
x0 = 1;
Fsumsquares = @(x)sum((specs.F_mean(x,xdata) - ydata).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[mean_param,ressquared,eflag,outputu] = fminunc(Fsumsquares,x0,opts)
% F_mean=fit(xdata',ydata','smoothingspline');
%%
ni_sd=ni_mean;%isnan(y_spike_sd) | y_spike_mean > 10*200;
cap_sd=y_spike_sd>quantile(y_spike_sd, 0.9);
xdata=x_current(~ni_sd & ~cap_sd);ydata=y_spike_sd(~ni_sd & ~cap_sd);

x0=[1 1];
Fsumsquares = @(x)sum((specs.F_sd(x,xdata) - ydata).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[sd_param,ressquared,eflag,outputu] =  fminunc(Fsumsquares,x0,opts);

% qbounds=quantile(ydata,[0.1 0.9]);
% sd_index= find(ydata > qbounds(1) & ydata < qbounds(2));
% sd_average=mean(ydata(sd_index));
% sd_average=median(ydata(sd_index));


% F_sd=fit(xdata',ydata','smoothingspline');
%% 
%% Obtain a grid of the mean and variance:

current_grid = specs.current_min:specs.current_gap:specs.current_max;
x_grid= current_grid;

spike_mean_grid=specs.F_mean(mean_param,x_grid);
spike_sd_grid=specs.F_sd(sd_param,x_grid);

% spike_mean_grid=F_mean(x_grid);
% spike_sd_grid=F_sd(x_grid);

% spike_sd_grid=sd_average.*ones(length(x_grid),1);

spike_curves=struct;
spike_curves.current=current_grid*specs.current_multiplier;
spike_curves.mean=spike_mean_grid;
spike_curves.sd=spike_sd_grid;

min_sd_index=find(spike_sd_grid<specs.sd_min);
spike_curves.sd(min_sd_index)=specs.sd_min;
max_sd_index=find(spike_sd_grid>specs.sd_max);
spike_curves.sd(max_sd_index)=specs.sd_max;

spike_curves.x_current=x_current;
spike_curves.y_spike_mean=y_spike_mean;
spike_curves.y_spike_sd=y_spike_sd;
