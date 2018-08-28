function [GP_samples] = GP_1D_tau_tuning(tau, locations,axis,pilot_data)

GP_params.mean_params=pilot_data.(axis).mean_params;
GP_params.var_params=pilot_data.(axis).var_params;
GP_params.tau=tau;
n_shapes=1;
GP_samples=draw_1D_GP(locations,n_shapes,GP_params);
