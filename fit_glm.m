function fit_params = fit_glm(stims_x, stims_t, spikes, init_vals)

% global g_stims_x
% global g_stims_t
% global g_spikes
% 
% g_stims_x = stims_x;
% g_stims_t = stims_t;
% g_spikes = spikes;

num_basis_funcs = 10;
num_spatial_pos = 121;

num_params = 1 + 2*num_basis_funcs + num_spatial_pos;

bs = basisFactory.makeSmoothTemporalBasis('raised cosine', .750, num_basis_funcs, @(x) 75);
basis = bs.B(1:end,:)';

% global g_basis
% g_basis = basis;

% x0= zeros(num_params,1) + 1e-9;
x0 = init_vals;
obj_fun = @(x) glm_log_likelihood(x, stims_x, stims_t, spikes, basis);
iterFun = @plot_params;

% setup.order = 2;
% setup.numvar = length(init_vals);
% setup.objective  = 'glm_log_likelihood';
% % setup.constraint = 'confun';

% adifuncs = adigatorGenFiles4Fmincon(setup);

options = optimoptions('fmincon','Algorithm','interior-point',...
                                'Display','iter','GradObj','on',...
                                'Diagnostics','on',...
                                'OutputFcn',iterFun); %
                            
% options = optimoptions('fminunc','Algorithm','trust-region', 'Display','iter','GradObj','on');

% options = optimset('Display','iter');

ub = Inf*ones(size(x0));                            
lb = -Inf*ones(size(x0));
lb(1:num_spatial_pos+num_basis_funcs+1) = 0;
% ub(1) = 1e-8;
ub(2:num_spatial_pos+1) = 1e-2;
                            

% fit_params = fminunc(obj_fun,x0,options)
fit_params = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
% fit_params = fminsearch(obj_fun,x0,options);



function stop = plot_params(x, optimValues, state)

stop = false;
if mod(optimValues.iteration,10) == 0
    fit_params = x;

    num_basis_funcs = 10;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', .750, num_basis_funcs, @(x) 75);
    basis = bs.B(1:end,:)';

    % build the filters from basis function
    stim_filter = zeros(1,size(basis,2));
    history_filter = zeros(1,size(basis,2));

    param_count = 1;
    baseline_rate = fit_params(param_count); param_count = param_count + 1;

    t = (1:1500)/1000;

    % build spatial layout of cell
    num_spatial_pos = 121;
    sptial_footprint = fit_params(param_count:param_count + num_spatial_pos-1);
    param_count = param_count + num_spatial_pos;

    % a_stim = fit_params(param_count); param_count = param_count + 1;
    % c_stim = fit_params(param_count); param_count = param_count + 1;

    for i = 1:num_basis_funcs
    %     phi = fit_params(param_count); param_count = param_count + 1;
        weight =  fit_params(param_count); param_count = param_count + 1;
        stim_filter = stim_filter + weight * basis(i,:); %raised_cosine(t,a_stim,c_stim,phi);
    end

    % a_hist = fit_params(param_count); param_count = param_count + 1;
    % c_hist = fit_params(param_count); param_count = param_count + 1;
    for i = 1:num_basis_funcs
    %     phi = fit_params(param_count); param_count = param_count + 1;
        weight =  fit_params(param_count); param_count = param_count + 1;
        history_filter = history_filter + weight * basis(i,:); %raised_cosine(t,a_hist,c_hist,phi);
    end

    figure(1)
    subplot(311)
    imagesc(reshape(sptial_footprint,11,11)')
    colorbar
    title(['Cell ' num2str(i)])

    subplot(312)
    plot(stim_filter)

    subplot(313)
    plot(history_filter)

    drawnow
end

function [c ceq] = confun(x)
c = 0; ceq = 0;