function [spike_curves, x_current, y_spike_mean] = get_spike_curves(curves_pilot_path,varargin)
% Revision notes:
%   1. This function should take a processed data from single-patch
%   experiments (that can include data from the shape model)
%   2. We will also need to remove the bad cell in the data, maybe create
%   another file for clean and procsssed single-patch data 

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

%% Read data from curves_pilot_path
load(curves_pilot_path); % An object named curves_pilot
x_current=curves_pilot.x_current;
y_spike_mean=curves_pilot.y_spike_mean;
y_spike_sd=curves_pilot.y_spike_sd;

if ~isempty(varargin) && ~isempty(varargin{1})
    specs= varargin{1};
else
    specs=struct;
    specs.time_max=200; % maximum spike time 
%     specs.F_mean = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3);
    %specs.F_mean = @(x,xdata) x(2)./(xdata + x(1)) + x(3); %min(y_spike_mean) + 
    % New mean function that reflects the saturation of the min spike time:
    specs.F_mean = @(x,xdata)  (x(2)./(xdata + x(1)) + x(3)).*(xdata<x(4))+(xdata>=x(4)).*(x(2)./(x(4) + x(1)) + x(3)) ; %min(y_spike_mean) + 
    
    specs.F_dev = @(x,xdata)  (x(2)./(xdata + x(1)).^2 + x(3)).*(xdata<x(4))+(xdata>=x(4)).*(x(2)./(x(4) + x(1)).^2 + x(3)) ; %min(y_spike_mean) + 
    specs.F_sd =  @(x,xdata)  (x(2)./(xdata + x(1)).^2 + x(3)).*(xdata<x(4))+(xdata>=x(4)).*(x(2)./(x(4) + x(1)).^2 + x(3)) ; %min(y_spike_mean) + 
    specs.current_multiplier=1e-3;
    specs.current_min=10;
    specs.current_max=10000;
    specs.current_gap=2;
    specs.sd_min=2;
    
    % The following will be rewritten with pilot data max
    specs.sd_max=15; 
    specs.dev_max=15; % bound the deviation as well
end

%% Call fmincon:
% y_spike_mean > specs.time_max 



ni_mean=isnan(y_spike_mean) |  x_current > specs.current_max | isnan(x_current);
assignin('base','y_spike_mean',y_spike_mean);
assignin('base','x_current',x_current);
% find(ni_mean)

% ni_mean=isnan(y_spike_mean) | x_current > 3500;
xdata=x_current(~ni_mean);ydata=y_spike_mean(~ni_mean);

% x0=[120 .01 min(y_spike_mean)];
x0 = [1 1 0 2000];
Fsumsquares = @(x)sum((specs.F_mean(x,xdata) - ydata).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
[mean_param,ressquared,eflag,outputu] = fminunc(Fsumsquares,x0,opts);
% F_mean=fit(xdata',ydata','smoothingspline');

%% Obtain a smooth curve the for errors in fitting a mean curve 
% Use a log transformation 
yfit=specs.F_mean(mean_param,xdata);
ydev=abs(yfit-ydata);
 
    specs.dev_max=max(ydev); 
ydev_log=log(ydev);
% index_dev= ydev > 1e3;

% x0=[120 .01 min(y_spike_mean)];
x0 = [1 -1 0 2000];
Fsumabs = @(x)sum((specs.F_dev(x,xdata) - ydev_log).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
[dev_param,ressquared,eflag,output] = fminunc(Fsumabs,x0,opts);
% F_mean=fit(xdata',ydata','smoothingspline');

%% Obtain a crude curve for SD (within condition variability)
% Use a log-transformation again
ni_sd=ni_mean | y_spike_sd<=0;%isnan(y_spike_sd) | y_spike_mean > 10*200;
% cap_sd=y_spike_sd>quantile(y_spike_sd, 0.9);
xdata=x_current(~ni_sd );
ysd=y_spike_sd(~ni_sd);
specs.sd_max=max(ysd); 
ysd_log=log(ysd);
x0 = [1 -1 0 2000];
Fsumsquares = @(x)sum((specs.F_sd(x,xdata) - ysd_log).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
[sd_param,ressquared,eflag,output] =  fminunc(Fsumsquares,x0,opts);

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
spike_sd_grid=exp(specs.F_sd(sd_param,x_grid));
% spike_dev_2_grid=specs.F_dev(dev_param,x_grid);
% index_neg=find(spike_dev_2_grid<0);
% spike_dev_2_grid(index_neg)=0;
spike_dev_grid=exp(specs.F_dev(dev_param,x_grid));

% spike_mean_grid=F_mean(x_grid);
% spike_sd_grid=F_sd(x_grid);

% spike_sd_grid=sd_average.*ones(length(x_grid),1);

spike_curves=struct;
spike_curves.current=current_grid*specs.current_multiplier;
spike_curves.mean=spike_mean_grid;
spike_curves.sd=spike_sd_grid;
spike_curves.dev=spike_dev_grid;

min_sd_index=find(spike_sd_grid<specs.sd_min);
spike_curves.sd(min_sd_index)=specs.sd_min;
max_sd_index=find(spike_sd_grid>specs.sd_max);
spike_curves.sd(max_sd_index)=specs.sd_max;

max_dev_index=find(spike_dev_grid>specs.dev_max);
spike_curves.dev(max_dev_index)=specs.dev_max;
min_dev_index=find(spike_dev_grid<0);
spike_curves.dev(min_dev_index)=0;

spike_curves.x_current=x_current;
spike_curves.y_spike_mean=y_spike_mean;
spike_curves.y_spike_sd=y_spike_sd;
spike_curves.mean_param=mean_param;
spike_curves.sd_param=sd_param;
spike_curves.dev_param=dev_param;
spike_curves.specs=specs;
 spike_curves.time_max=specs.time_max;
% %%
% scatter(x_current,y_spike_mean)
% hold on;
% plot(x_grid,spike_mean_grid) 
% 
% ylim([min(y_spike_mean)-20 max(y_spike_mean)+20])
% %%
% scatter(x_current,y_spike_sd)
% hold on;
% plot(x_grid,spike_sd_grid) 
% ylim([min(y_spike_sd)-20 max(y_spike_sd)+20])
% hold on;
% %%
% scatter(xdata,ydev,'MarkerFaceColor','red');
% hold on;
% plot(x_grid,spike_dev_grid) 
% ylim([min(ydev) max(ydev)])
