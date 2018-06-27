function  [mean_params, var_params,neurons]=pre_processing(neurons,params)
%%
n_cell = length(neurons);
for i_cell = 1:n_cell
    [ugrid,~,ua]=unique(neurons(i_cell).stim_grid);
    neurons(i_cell).unique_grid = ugrid;
    neurons(i_cell).avg_current=zeros(length(ugrid),1);
    neurons(i_cell).raw_sq_deviation=zeros(length(ugrid),1);
    neurons(i_cell).mean_current=zeros(length(neurons(i_cell).raw_current),1);
    neurons(i_cell).sq_deviation=zeros(length(ugrid),1);
    neurons(i_cell).mean_power=zeros(length(ugrid),1);% not necessary since power is the same per loc & cell
    
    if params.power_scaling
        neurons(i_cell).current=...
            neurons(i_cell).raw_current./(neurons(i_cell).power);
    else
        neurons(i_cell).current=...
            neurons(i_cell).raw_current;
    end
    for i_grid = 1:length(ugrid)
        indices = find(ua==i_grid);
        % these are defined on the unique grid point
        neurons(i_cell).avg_current(i_grid) =mean(neurons(i_cell).current(indices));
        neurons(i_cell).mean_power(i_grid) =mean(neurons(i_cell).power(indices));
        neurons(i_cell).raw_sq_deviation(i_grid) =var(neurons(i_cell).raw_current(indices));
        % the following is defined for each data point:
        neurons(i_cell).mean_raw_current(indices)=mean(neurons(i_cell).raw_current(indices));
    end
    max_current = max(neurons(i_cell).avg_current);
    neurons(i_cell).scaled_current = neurons(i_cell).current/max_current;
    neurons(i_cell).scale=max_current;
    
    % Assign a sigma for each data point:
    neurons(i_cell).noise_sigma =  ones(length(neurons(i_cell).current),1)*sqrt(mean(neurons(i_cell).sq_deviation))/(neurons(i_cell).scale);
end

% Sigma ~ current model
if params.sigma_current_model
    noise_var=[];mean_current =[];
    
    for i_cell =1:n_cell
        noise_var = [noise_var neurons(i_cell).raw_sq_deviation'];
        mean_current = [mean_current neurons(i_cell).avg_current'];
    end
    dims=[size(noise_var,1)*size(noise_var,2) 1];
    Y=reshape(sqrt(noise_var), dims);
    X= [ones(dims) reshape(mean_current, dims)];
%     X= [reshape(mean_current, dims)];

% 
fun=@(x) sum(( Y-X*x').^2);
    A = [];
b = [];
Aeq = [];
beq = [];
lb = [params.sd_min 0];
ub = [Inf Inf];
x0 = [0,0];
    [est,fval,exitflag,output] =fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
    
    for i_cell = 1:n_cell
        % NOTE: the current code doest not include scaling by power!
        scaling_factor= ones(length(neurons(i_cell).current),1)*neurons(i_cell).scale;
        if params.power_scaling
            scaling_factor= neurons(i_cell).power.*scaling_factor;
        end
        neurons(i_cell).noise_sigma =  slope*neurons(i_cell).mean_raw_current'./scaling_factor;
        neurons(i_cell).slope=est(2);
        neurons(i_cell).intercept=est(1);
        
    end
end
%% Find centers using isometic regresson:
% plot(neurons(1).stim_grid,neurons(1).scaled_current)
% 1 to ngrid +1
for i_cell = 1:n_cell
    [unique_stim,unique_ia,unique_ic]=unique(neurons(i_cell).stim_grid);
    mean_current =zeros(1,max(unique_ic));
    for i = 1:max(unique_ic)
        mean_current(i)= mean(neurons(i_cell).scaled_current(unique_ic==i));
    end
    
    n_grid = length( mean_current);
    [left]=prefixiso(n_grid,mean_current);
    [right]=prefixiso(n_grid,flip(mean_current));
    errorl=left.error(2:end);
    errorr=flip(right.error(2:end));
    totalerror=zeros(n_grid,1);
    for i=2:n_grid
        totalerror(i) = errorl(i-1)+errorr(i);
    end
    
    totalerror(1)=max(totalerror);
    [~,I_center]=min(totalerror);
    neurons(i_cell).initial_shift = unique_stim(I_center);
end

for i_cell = 1:n_cell
    neurons(i_cell).adjusted_grid = neurons(i_cell).stim_grid-neurons(i_cell).initial_shift;
end
%% Fit the mean funtion from the centered data:




X=[neurons(:).adjusted_grid]';
Y=[neurons(:).scaled_current]';

if params.symmetric
    X=[X; -X];
    Y=[Y; Y];
end


params.X=X;
params.Y=Y;

fixed_grid = -(params.boundary+params.buffer ):0.1:(params.boundary+params.buffer );
Xstar=fixed_grid';


[pred_mean, post_mean,prior_var] = get_GP_boundary(Xstar, params);
mean_params.grid=fixed_grid;
mean_params.values=pred_mean;
mean_params.prior_var=prior_var;
mean_params.data.Y = params.Y;
mean_params.data.X = params.X;

%     figure(1)
%     scatter(mean_params.data.X,post_mean)
%     hold on;
%     plot(mean_params.grid,mean_params.values)
%     hold on;
%     scatter(mean_params.data.X,mean_params.data.Y)
%% Now estimate the variance over the full shape space:
observed_values =  params.Y;
fitted_values = post_mean;
sqdev=(observed_values-fitted_values).^2;

% change the response variable to squared deviation
params.Y=sqdev;

[pred_var, post_var,prior_var] = get_GP_boundary(Xstar, params);
var_params.grid=fixed_grid;
var_params.values=pred_var;
var_params.prior_var=prior_var;

var_params.data.Y = params.Y;
var_params.data.X = params.X;

%%
for i_cell = 1:n_cell
    neurons(i_cell).fit_gain=params.fit_gain;
end


