function  [pilot_data]=pre_processing_shape(pilot_data,params)
%% Reformat the data and learn the shapes 
% 1) find shift
% Based on the estimated shifts
% 2) estimate mean function
% 3) estimate variance function
zero_bound = 1e-10;
axis_list = fieldnames(pilot_data);
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    neurons=pilot_data.(ax).neurons;
    n_cell = length(neurons);
    for i_cell = 1:n_cell
        if strcmp(ax,'xy')
            [ugrid,~,ua]=unique(neurons(i_cell).stim_grid','rows');
        else
            [ugrid,~,ua]=unique(neurons(i_cell).stim_grid);
        end
        neurons(i_cell).unique_grid = ugrid;
        neurons(i_cell).avg_current=zeros(length(ugrid),1);
        neurons(i_cell).avg_raw_current=zeros(length(ugrid),1);
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
            neurons(i_cell).avg_raw_current(i_grid) =mean(neurons(i_cell).raw_current(indices));
            
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
    pilot_data.(ax).neurons=neurons;
end
%% Set the nugget:
pilot_data.nugget=params.nugget;
%% Learn the center shifts 
% Use the heuristic approach by picking the peak and the location of nuclei
if params.find_shifts
 % 1 to ngrid +1
    for i_ax = 1:length(axis_list)
        ax=axis_list{i_ax};
        neurons=pilot_data.(ax).neurons;
        n_cell = length(neurons);
        if ~strcmp(ax,'xy')
            for i_cell = 1:n_cell
                [unique_stim,unique_ia,unique_ic]=unique(neurons(i_cell).stim_grid);
                mean_current =zeros(1,max(unique_ic));
                for i = 1:max(unique_ic)
                    mean_current(i)= mean(neurons(i_cell).scaled_current(unique_ic==i));
                end
                [~, I_center]=max(mean_current);
                neurons(i_cell).initial_shift = unique_stim(I_center);
            end
            for i_cell = 1:n_cell
                neurons(i_cell).adjusted_grid = neurons(i_cell).stim_grid-neurons(i_cell).initial_shift;
                neurons(i_cell).find_shift=params.find_shifts;
            end
        pilot_data.(ax).neurons=neurons;
        else
            for i_cell = 1:n_cell
                [unique_stim,unique_ia,unique_ic]=unique(neurons(i_cell).stim_grid', 'rows');
                mean_current =zeros(1,max(unique_ic));
                for i = 1:max(unique_ic)
                    mean_current(i)= mean(neurons(i_cell).scaled_current(unique_ic==i));
                end
                [~, I_center]=max(mean_current);
                neurons(i_cell).initial_shift = unique_stim(I_center,:);
            end
               for i_cell = 1:n_cell
                neurons(i_cell).adjusted_grid = neurons(i_cell).stim_grid-transpose(neurons(i_cell).initial_shift)*ones(1,size(neurons(i_cell).stim_grid,2));
                neurons(i_cell).find_shift=params.find_shifts;
            end
        pilot_data.(ax).neurons=neurons;
        end
    end
   else % Take 0 as the mode
    for i_ax = 1:length(axis_list)
        ax=axis_list{i_ax};
        neurons=pilot_data.(ax).neurons;
        n_cell = length(neurons);
        for i_cell = 1:n_cell
            neurons(i_cell).initial_shift=0;
            neurons(i_cell).adjusted_grid = neurons(i_cell).stim_grid;
            neurons(i_cell).find_shift=params.find_shifts;
        end
        pilot_data.(ax).neurons=neurons;
    end
end
%% Fit the mean funtion from the centered data:
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    neurons=pilot_data.(ax).neurons;
    if ~strcmp(ax,'xy')
        X=[neurons(:).adjusted_grid]';
        Y=max(zero_bound,[neurons(:).scaled_current]);
        Y= params.inv_excite(Y)';
        if params.symmetric
            X=[X; -X];
            Y=[Y; Y];
        end
        params.(ax).X=X;
        params.(ax).Y=Y;
        params.(ax).prior_mean=params.inv_excite(zero_bound);
        
        fixed_grid = -(params.(ax).boundary+params.(ax).buffer ):0.1:(params.(ax).boundary+params.(ax).buffer );
        Xstar=fixed_grid';
        
        [pred_mean, post_mean,prior_var] = get_GP_boundary(Xstar, params.(ax));
%      post_mean = smooth(x,y,0.1,'loess');
        mean_params.grid=fixed_grid;
        mean_params.values=pred_mean;
        mean_params.prior_var=prior_var;
        mean_params.data.Y = params.(ax).Y;
        mean_params.data.fitted=post_mean;
        mean_params.data.X = params.(ax).X;
        mean_params.inv_excite=  params.inv_excite;  mean_params.excite=  params.excite;
        pilot_data.(ax).mean_params=mean_params;
        
        neurons=pilot_data.(ax).neurons;
        observed_values =  params.(ax).Y;
        fitted_values =mean_params.data.fitted;
        sqdev=(observed_values-fitted_values).^2;

        % change the response variable to squared deviation
        params.(ax).Y=sqdev;
        if isfield(params.(ax),'distance_factor_var')
            params.(ax).distance_factor=params.(ax).distance_factor_var;
        end
        if isfield(params.(ax),'tau_var')
            params.(ax).tau=params.(ax).tau_var;
        end
          params.(ax).prior_mean=zero_bound;
      
        [pred_var, post_var,prior_var] = get_GP_boundary(Xstar, params.(ax));
        var_params.grid=fixed_grid;
        var_params.values=max(zero_bound,pred_var);
        var_params.prior_var=prior_var;
        
        var_params.data.Y = params.(ax).Y;
        var_params.data.X = params.(ax).X;
        var_params.inv_excite=  params.inv_excite;  var_params.excite=  params.excite;
        
        pilot_data.(ax).var_params=var_params;
        pilot_data.(ax).tau=params.(ax).tau;
    else

        X=[neurons(:).adjusted_grid]';
        X=X(:,1:2);
        Y=max(zero_bound,[neurons(:).scaled_current]);
        Y=params.inv_excite(Y)';
        params.(ax).X=X;
        params.(ax).Y=Y;
        params.(ax).prior_mean=params.inv_excite(zero_bound);
      
          
        fixed_grid_x = -(params.(ax).x.boundary+params.(ax).x.buffer ):3:(params.(ax).x.boundary+params.(ax).x.buffer );
        fixed_grid_y = -(params.(ax).y.boundary+params.(ax).y.buffer ):3:(params.(ax).y.boundary+params.(ax).y.buffer );
        [temp1, temp2]=meshgrid(fixed_grid_x,fixed_grid_y);
        Xstar=[reshape(temp1,[],1)  reshape(temp2,[],1)];
        %     Xstar=neurons(1).unique_grid;
        [pred_mean, post_mean,prior_var] = get_2DGP_boundary(Xstar, params.(ax));
        mean_params.grid=Xstar;
        mean_params.values=pred_mean;
        mean_params.prior_var=prior_var;
        mean_params.data.Y = params.(ax).Y;
        mean_params.data.fitted=post_mean;
        mean_params.data.X = params.(ax).X;
        mean_params.inv_excite=  params.inv_excite;  mean_params.excite=  params.excite;
        pilot_data.(ax).mean_params=mean_params;
        
        
        neurons=pilot_data.(ax).neurons;
        
        observed_values =  params.(ax).Y;
        fitted_values = pilot_data.(ax).mean_params.data.fitted;
        sqdev=(observed_values-fitted_values).^2;
        % change the response variable to squared deviation
        params.(ax).Y=sqdev;
          params.(ax).prior_mean=zero_bound;
      
        [pred_var, post_var,prior_var] = get_2DGP_boundary(Xstar, params.(ax));
        var_params.grid=Xstar;
        var_params.values=max(zero_bound,pred_var);
        var_params.prior_var=prior_var;
        var_params.data.Y = params.(ax).Y;
        var_params.data.X = params.(ax).X;
        var_params.inv_excite=  params.inv_excite;  var_params.excite=  params.excite;
        
        pilot_data.(ax).var_params=var_params;
        
    pilot_data.(ax).tau=[params.(ax).x.tau params.(ax).y.tau];
    end
end
%     figure(1)
%     scatter(mean_params.data.X,post_mean)
%     hold on;
%     plot(mean_params.grid,mean_params.values)
%     hold on;
%     scatter(mean_params.data.X,mean_params.data.Y)
%% Update the neurons for each dimension 
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    neurons=pilot_data.(ax).neurons;
    n_cell = length(neurons);
    for i_cell = 1:n_cell
        neurons(i_cell).fit_gain=params.fit_gain;
    end
    pilot_data.(ax).neurons=neurons;
end

%% Summarize the mean and variance of shifts 
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    pilot_data.(ax).shift_params =struct;
    if ~strcmp(ax,'xy')
        
        pilot_data.(ax).shift_params.data = [pilot_data.(ax).neurons(:).initial_shift];
        pilot_data.(ax).shift_params.mean = mean([pilot_data.(ax).neurons(:).initial_shift]);
        pilot_data.(ax).shift_params.var = var([pilot_data.(ax).neurons(:).initial_shift]);
    else
        tmp=reshape([pilot_data.xy.neurons.initial_shift],3,[]);
        pilot_data.(ax).shift_params.data= tmp(1:2,:);
        pilot_data.(ax).shift_params.mean= mean(tmp(1:2,:)')';
        pilot_data.(ax).shift_params.var= cov(tmp(1:2,:)')';
    end
end
%% Convolute the GP mean functions using the distributions of shifts 
% i.e., marginalizing the shift to allow for greater robustness 
% moved to this script since it took a while to run 
axis_list = fieldnames(pilot_data);
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    mean_params=  pilot_data.(ax).mean_params;
    shift_params=  pilot_data.(ax).shift_params;
    tmp=mean_params.values;
         if ~strcmp(ax,'xy')
            shift_grid = -3:0.5:3; % This range should contain almost all probability mass for a standard normal r.v.
            prob_weight=normpdf(shift_grid);prob_weight=prob_weight/sum(prob_weight); prob_weight=prob_weight';
            shift_grid = shift_grid*sqrt(shift_params.var)+shift_params.mean; % Transform it given the shift distribution 
         else
           [X Y]=meshgrid(-3:0.5:3,-3:0.5:3);
           shift_grid  = [ reshape(X,[],1) reshape(Y,[],1)];
            prob_weight=normpdf(shift_grid);prob_weight=prob_weight/sum(prob_weight);
           shift_grid=shift_grid*sqrtm(shift_params.var) + ones(size(shift_grid,1),1)*(shift_params.mean');
         end
    for i_grid = 1:length(mean_params.values)
            if ~strcmp(ax,'xy')
             X_r = mean_params.grid(i_grid)+shift_grid;
             X_r=X_r';
            else
                X_r = mean_params.grid(i_grid,:)+shift_grid;
            end
        tmp(i_grid)=sum(quick_match(X_r,mean_params).*prob_weight);
        
    end
     pilot_data.(ax).mean_params.convoluted_values=tmp;
end

%% Gather data for learning the sigma - current relationships 
% We can acutally ignore this property for building the prior!
if params.sigma_current_model
    if params.joint_sigma % fit one single model across x,y, and z
        noise_var=[];mean_current =[];
        for i_ax = 1:length(axis_list)
            ax=axis_list{i_ax};
            neurons=pilot_data.(ax).neurons;
            n_cell = length(neurons);
            for i_cell =1:n_cell
                noise_var = [noise_var neurons(i_cell).raw_sq_deviation'];
                mean_current = [mean_current neurons(i_cell).avg_raw_current'];
            end
        end
        dims=[size(noise_var,1)*size(noise_var,2) 1];
        Y=reshape(sqrt(noise_var), dims);
        X= [ones(dims) reshape(mean_current, dims)];
        %     X= [reshape(mean_current, dims)];
        fun=@(x) sum(( Y-X*x').^2);
        A = [];b = [];Aeq = [];beq = []; lb = [params.sd_min 0];
        ub = [Inf Inf]; x0 = [0,0];
        [est,fval,exitflag,output] =fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
        for i_ax = 1:length(axis_list)
            ax=axis_list{i_ax};
            neurons=pilot_data.(ax).neurons;
            n_cell = length(neurons);
            for i_cell = 1:n_cell
                % NOTE: the current code doest not include scaling by power!
                scaling_factor= ones(1,length(neurons(i_cell).current))*neurons(i_cell).scale;
                if params.power_scaling
                    scaling_factor= neurons(i_cell).power.*scaling_factor;
                end
                neurons(i_cell).noise_sigma = (est(1)+est(2)*neurons(i_cell).mean_raw_current)./scaling_factor;
                neurons(i_cell).slope=est(2);
                neurons(i_cell).intercept=est(1);
            end
            pilot_data.(ax).neurons=neurons;
        end
        pilot_data.sd_current =struct;
        pilot_data.sd_current.X=X;
        pilot_data.sd_current.Y=Y;
        pilot_data.sd_current.beta=est;
        
    else % fit the sigma-current model separately
        % Sigma ~ current model (use all data)
        for i_ax = 1:length(axis_list)
            ax=axis_list{i_ax};
            neurons=pilot_data.(ax).neurons;
            n_cell = length(neurons);
            noise_var=[];mean_current =[];
            for i_cell =1:n_cell
                noise_var = [noise_var neurons(i_cell).raw_sq_deviation'];
                mean_current = [mean_current neurons(i_cell).avg_raw_current'];
            end
            dims=[size(noise_var,1)*size(noise_var,2) 1];
            Y=reshape(sqrt(noise_var), dims);
            X= [ones(dims) reshape(mean_current, dims)];
            %     X= [reshape(mean_current, dims)];
            fun=@(x) sum(( Y-X*x').^2);
            A = []; b = [];Aeq = [];
            beq = []; lb = [params.sd_min 0];
            ub = [Inf Inf]; x0 = [0,0];
            [est,fval,exitflag,output] =fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
            for i_cell = 1:n_cell
                % NOTE: the current code doest not include scaling by power!
                scaling_factor= ones(1,length(neurons(i_cell).current))*neurons(i_cell).scale;
                if params.power_scaling
                    scaling_factor= neurons(i_cell).power.*scaling_factor;
                end
                neurons(i_cell).noise_sigma = (est(1)+est(2)*neurons(i_cell).mean_raw_current)./scaling_factor;
                neurons(i_cell).slope=est(2);
                neurons(i_cell).intercept=est(1);
            end
            pilot_data.(ax).neurons=neurons;
        end
    end
end


%% Find centers using isometic regresson:
% Not in use 
% if params.find_shifts_isometric % Fit isometric regression to find the mode.
%     plot(neurons(1).stim_grid,neurons(1).scaled_current)
%     1 to ngrid +1
%     for i_ax = 1:length(axis_list)
%         ax=axis_list{i_ax};
%         neurons=pilot_data.(ax).neurons;
%         n_cell = length(neurons);
%         
%         for i_cell = 1:n_cell
%             [unique_stim,unique_ia,unique_ic]=unique(neurons(i_cell).stim_grid);
%             mean_current =zeros(1,max(unique_ic));
%             for i = 1:max(unique_ic)
%                 mean_current(i)= mean(neurons(i_cell).scaled_current(unique_ic==i));
%             end
%             n_grid = length( mean_current);
%             [left]=prefixiso(n_grid,mean_current);
%             [right]=prefixiso(n_grid,flip(mean_current));
%             errorl=left.error(2:end);
%             errorr=flip(right.error(2:end));
%             totalerror=zeros(n_grid,1);
%             for i=2:n_grid
%                 totalerror(i) = errorl(i-1)+errorr(i);
%             end
%             totalerror(1)=max(totalerror);
%             [~,I_center]=min(totalerror);
%             neurons(i_cell).initial_shift = unique_stim(I_center);
%         end
%         for i_cell = 1:n_cell
%             neurons(i_cell).adjusted_grid = neurons(i_cell).stim_grid-neurons(i_cell).initial_shift;
%             neurons(i_cell).find_shift=params.find_shifts;
%         end
%         pilot_data.(ax).neurons=neurons;
%     end
% else % Take 0 as the mode
%     for i_ax = 1:length(axis_list)
%         ax=axis_list{i_ax};
%         neurons=pilot_data.(ax).neurons;
%         n_cell = length(neurons);
%         for i_cell = 1:n_cell
%             neurons(i_cell).initial_shift=0;
%             neurons(i_cell).adjusted_grid = neurons(i_cell).stim_grid;
%             neurons(i_cell).find_shift=params.find_shifts;
%         end
%         pilot_data.(ax).neurons=neurons;
%     end
% end