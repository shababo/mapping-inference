
%% An alternative way to generate the samples 
% The assignments depend on the gamma and f
% Last update: Jan. 29th 2017. 
% Note: The variational distribution of $a$ is problematic. 


% This means that we might want to merge the likelihood estimation and
% sampling stage
%% Specify the prior distribution 
prior_parameters = struct('prob',0.2,'mu0',5,'sigma0',1);
inhibit_period = 1;
rho = 0.1;
%% Initialization
prob_initial = 0; %logit(1)
mean_initial = 6;
sd_initial = 1;

var_params = struct('prob',[],'mean',[],'sd',[],...
    'phi',[]);
% Synaptic parameters: p, m, s
% Assignments parameters: phi 

var_params.prob = prob_initial.*ones(p,1);
var_params.mean = mean_initial.*ones(p,1);
var_params.sd = sd_initial.*ones(p,1);

% Initialize phi:
% We need to create a p-vector of zero for each event 
% Note that phi_0 is set to be 0 for identifiability

var_params.phi = cell(N,1);

for ni = 1:N
    if size(postsynaptic_events{ni}.event_times,2) > 0
        var_params.phi{ni} = zeros(p,size(postsynaptic_events{ni}.event_times,2));
    end
end

% Number of samples:
L = 50; %might need a lot of samples

%% Start of iteration:

%% Draw samples from q()
% We sequentially draw
% 1, gamma, mu
% 2, assignments at increasing time

% Grid for estimating the intensity
MCsamples = cell([L 1]);
for l = 1:L
    MCsamples{l} = struct('gamma',[],'mu',[],'assignments',[]);
    % Draw gamma
    MCsamples{l}.gamma = unifrnd(0,1,p,1) < exp(var_params.prob)./(1+exp(var_params.prob ));
    % Draw mu
    MCsamples{l}.mu = normrnd(var_params.mean, var_params.sd);
    % Note that the exact value of mu does not matter if gamma = 0
    % To avoid errors, we set mu to be zero if gamma = 0
    MCsamples{l}.mu(MCsamples{l}.gamma==0) = 0;
    
    % Draw assignments
    
    MCsamples{l}.assignments = cell([N 1]);
    MCsamples{l}.assignments_intensity = cell([N 1]);
    
    
    for ni = 1:N
        events_this_trial = postsynaptic_events{ni}.event_times;
        MCsamples{l}.assignments{ni} = events_this_trial;
        event_counts = size(events_this_trial,2);
        
        if event_counts > 0
            for event_i = 1:event_counts
                
                % Evaluate the current intensity (also their previous values:
                nonzero_intensities = zeros(p+1,1);
                nonzero_intensities(1) = spontaneous_intensity;% the intensity for spontaneous event is always nonzero
                
                for cell_i = 1:p
                    if parameters_tmp.gamma(cell_i) == 0 % if this cell is not connected
                        nonzero_intensities(cell_i+1)=0;
                    else
                        % Evaluate the intensity function up to this event time
                        
                        % TO-DO:
                        % Create a self-history object
                        % Evaluate the function up to that point..
                        f_current = 0;
                        if event_i == 1
                            f_current =  induced_intensity{ni}(cell_i);
                        else
                            events_idx = find(MCsamples{l}.assignments{ni} == cell_i);
                            if length(events_idx)>0
                                if max(events_this_trial(events_idx))+inhibit_period - events_this_trial(event_i) < 0
                                    f_current =  induced_intensity{ni}(cell_i);
                                end
                            end
                        end
                        nonzero_intensities(cell_i+1)= f_current;
                    end
                end
                prob_tmp = [1; exp(var_params.phi{ni}(:,event_i))];
                prob_tmp = prob_tmp .* (nonzero_intensities > 0);
                prob_tmp = transpose(prob_tmp./sum(prob_tmp));
                assign_tmp =mnrnd(1,prob_tmp);
                MCsamples{l}.assignments{ni}(event_i) =  find(assign_tmp)-1;
                MCsamples{l}.assignments_intensity{ni}(event_i) = nonzero_intensities(find(assign_tmp));
            end
        end
    end
end
%% Evaluate the stochastic gradient
Ngrid = 101;

fgrid = ([1:Ngrid] - 1).*Tmax./(Ngrid-1);
ftmp = zeros(Ngrid,1);
% We only use this grid to evaluate the integral
% The actuall intensity at each event can be very sensitive in time
% So we estimate them when drawing the assignments 

for l = 1:L
    
    parameters_tmp = MCsamples{l};
    % Calculate the log joint distribution and log q for each sample
    
    % Calculate the (log) likelihood for each trial
    loglklh = zeros(N,1);
    loglklh_p = zeros(p,1);
    for ni = 1:N
        % Evaluate the background events
        
        loglklh(ni) = - spontaneous_intensity*Tmax;
        spontaneous_mu = 6;
        events_idx = find(MCsamples{l}.assignments{ni} == 0);
        
        if length(events_idx)>0
            event_times_tmp = postsynaptic_events{ni}.event_times(events_idx);
            event_size_tmp = postsynaptic_events{ni}.event_size(events_idx);
            size_weights = normpdf(event_size_tmp,spontaneous_mu,size_sigma);
            loglklh(ni) = loglklh(ni)+sum(log(size_weights)) + length(events_idx)*log(spontaneous_intensity);
        end
        
        for cell_i = 1:p
            if parameters_tmp.gamma(cell_i) == 1
                % Extrat the assigned event for this cell
                events_idx = find(MCsamples{l}.assignments{ni} == cell_i);
                
                if length(events_idx)>0
                    event_times_tmp = postsynaptic_events{ni}.event_times(events_idx);
                    event_size_tmp = postsynaptic_events{ni}.event_size(events_idx);
                    
                    % Call a function here
                    for igrid = 1:Ngrid
                        most_recent_event = max( event_times_tmp(event_times_tmp <fgrid(igrid)));
                        
                        if length(most_recent_event)>0
                            if fgrid(igrid) - most_recent_event < 1 % 1 is a parameter
                                ftmp(igrid) = 0;
                            else
                                ftmp(igrid) = induced_intensity{ni}(cell_i);
                            end
                        else
                            ftmp(igrid) = induced_intensity{ni}(cell_i);
                        end
                    end
                    
                    % Integrate the intensity
                    loglklh_p(cell_i) = - sum(ftmp)*Tmax/Ngrid;
                    % Evaluate the density of the assigned events
                    size_weights = normpdf(event_size_tmp,parameters_tmp.mu(cell_i),size_sigma);
                    for event_i = 1:length(events_idx)
                        loglklh_p(cell_i) = loglklh_p(cell_i)+ log(ftmp(max(find(fgrid < event_times_tmp(event_i))))) ...
                            + log(size_weights(event_i));
                    end
                    
                else % No events are assigned to this cell
                    loglklh_p(cell_i) = - induced_intensity{ni}(cell_i)*Tmax;
                end
                
            else
                loglklh_p(cell_i) =0;
            end
            
        end
        loglklh(ni) = loglklh(ni) + sum(loglklh_p);
    end
    
    % Calculate the (log) prior distribution
    logprior  =0;
    for cell_i = 1:p
        if parameters_tmp.gamma(cell_i) == 1
            
            logprior  = logprior + log(prior_parameters.prob) + log(normpdf(parameters_tmp.mu(cell_i),...
                prior_parameters.mu0,prior_parameters.sigma0));
        else
            logprior  = logprior + log(1-prior_parameters.prob);
        end
        
    end
    
    % Evaluate log q()
    logq = 0;
    % log q(mu,gamma )
    for cell_i = 1:p
        if parameters_tmp.gamma(cell_i) == 1
            
            logq = logq + log(exp(var_params.prob(cell_i))/(1+exp(var_params.prob(cell_i))))+...
                log(normpdf(parameters_tmp.mu(cell_i), ...
                var_params.mean(cell_i), var_params.sd(cell_i)));
        else
            logq = logq + log(1/(1+exp(var_params.prob(cell_i))));
        end
    end
    
    % log q(a| gamma mu)
    logq_a = 0;
    for ni = 1:N
        if length(parameters_tmp.assignments{ni}) > 0
            for event_i = 1:length(parameters_tmp.assignments{ni})
                asgn_tmp = [1; exp(var_params.phi{ni}(:,event_i))];
                logq_a =logq_a+ log(asgn_tmp(parameters_tmp.assignments{ni}(event_i)+1)/sum(asgn_tmp));
            end
        end
        
    end
    logq = logq + logq_a;
    
    % Divergence
    total_divergence = sum(loglklh)+logprior - logq;
    
    % Evaluate the gradient w.r.t. each parameter
    gradient_wrt = var_params;
    
    
    for cell_i = 1:p
        gradient_wrt.prob(cell_i) = parameters_tmp.gamma(cell_i) - ...
            exp(var_params.prob(cell_i))/ (1+exp(var_params.prob(cell_i)));
        gradient_wrt.prob(cell_i) = gradient_wrt.prob(cell_i)*total_divergence;
        if parameters_tmp.gamma(cell_i)  == 1
            gradient_wrt.mean(cell_i) = - (var_params.mean(cell_i)-parameters_tmp.mu(cell_i))*...
            total_divergence/var_params.sd(cell_i)^2;
            
            gradient_wrt.sd(cell_i) = -1/var_params.sd(cell_i)+(var_params.mean(cell_i)-parameters_tmp.mu(cell_i))^2/...
                var_params.sd(cell_i)^3*total_divergence;
        else
            gradient_wrt.mean(cell_i)= 0;
            gradient_wrt.sd(cell_i)=0;
        end
        
    end
    for ni = 1:N
        if size(postsynaptic_events{ni}.event_times,2) > 0
            for event_i = 1:size(postsynaptic_events{ni}.event_times,2)
                gradient_wrt.phi{ni}(:,event_i) = - exp(var_params.phi{ni}(:,event_i))/(1+sum(var_params.phi{ni}(:,event_i)));
                if parameters_tmp.assignments{ni}(event_i) ==0
                    % do nothing
                else
                    idx_tmp =parameters_tmp.assignments{ni}(event_i);
                    gradient_wrt.phi{ni}(idx_tmp,event_i)  = 1+gradient_wrt.phi{ni}(idx_tmp,event_i);
                end
                gradient_wrt.phi{ni}(:,event_i) = gradient_wrt.phi{ni}(:,event_i)*total_divergence;
            end
        end
    end
    
    % Merge the results 
    if l == 1
        mean_gradient_wrt.mean = gradient_wrt.mean;
        mean_gradient_wrt.sd = gradient_wrt.sd;
        mean_gradient_wrt.phi =  gradient_wrt.phi;
    elseif l== L
        mean_gradient_wrt.mean = mean_gradient_wrt.mean+gradient_wrt.mean;
        mean_gradient_wrt.sd = mean_gradient_wrt.sd+gradient_wrt.sd;
        for ni = 1:N
            mean_gradient_wrt.phi{ni} = mean_gradient_wrt.phi{ni}+ gradient_wrt.phi{ni};
        end
    else
        mean_gradient_wrt.mean = rho*(mean_gradient_wrt.mean+gradient_wrt.mean)/L;
        mean_gradient_wrt.sd = rho*(mean_gradient_wrt.sd+gradient_wrt.sd)/L;
        for ni = 1:N
            mean_gradient_wrt.phi{ni} = rho*(mean_gradient_wrt.phi{ni}+ gradient_wrt.phi{ni})/L;
        end
    end
    l
    mean_gradient_wrt.mean
end
