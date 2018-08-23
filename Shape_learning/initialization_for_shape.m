 function [variational_params, prior_params,inference_params]=initialization_for_shape(i_ax,data_this_axis)
    ini_GP_tau=log(unifrnd(9,200));
    if i_ax == 3
        ini_GP_tau=log(unifrnd(25,500));
    end
    
    shift_sigma=2;
    tau_sigma=2;
    gain_sigma=2;
    ini_gain = log(unifrnd(1,20));
    
    neurons=data_this_axis.neurons;
    mean_params=data_this_axis.mean_params;
    var_params=data_this_axis.var_params;
    n_cell=length(neurons);
    
    clear('variational_params')
    variational_params(n_cell)=struct;
    for i_cell = 1:n_cell
        if neurons(1).find_shift
        variational_params(i_cell).shift.dist='normal';
        variational_params(i_cell).shift.type='individual';
        variational_params(i_cell).shift.mean=0;
        variational_params(i_cell).shift.log_sigma=log(shift_sigma);
        end
        
        variational_params(i_cell).GP_tau.dist='log-normal';
        variational_params(i_cell).GP_tau.type='common';
        variational_params(i_cell).GP_tau.mean=ini_GP_tau;
        variational_params(i_cell).GP_tau.log_sigma=log(tau_sigma);
        
        if neurons(1).fit_gain
            variational_params(i_cell).gain.dist='log-normal';
            variational_params(i_cell).gain.type='individual';
            variational_params(i_cell).gain.mean=ini_gain;
            variational_params(i_cell).gain.log_sigma=log(gain_sigma);
        end
    end
    prior_params=variational_params;
   
    inference_params=struct;
    inference_params.convergence_threshold=1e-3;
    inference_params.MCsamples_for_gradient=50;
    inference_params.step_size=1;
    inference_params.step_size_max=1;
    inference_params.maxit=2500;
    inference_params.eig_epsilon=1e-6;
    inference_params.mean_func=struct;
        inference_params.mean_func.func=@quick_match;
        inference_params.mean_func.params=mean_params;
    inference_params.var_func=struct;
        inference_params.var_func.func=@quick_match;
        inference_params.var_func.params=var_params;
    inference_params.prior_opt=false;
    inference_params.fit_gain=neurons(1).fit_gain;
    
        