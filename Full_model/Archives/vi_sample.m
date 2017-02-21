%% Apply variational inference on the toy data set

% Initialize variational parameters:
prob_initial = -1.386; %logit(0.2)
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
MCsamples = cell([L 1]);
for l = 1:L
   MCsamples{l} = struct('gamma',[],'mu',[],'assignments',[]);
   % Draw gamma
   MCsamples{l}.gamma = unifrnd(0,1,p,1) < exp(var_params.prob)./(1+exp(var_params.prob ));
   % Draw mu
   MCsamples{l}.mu = normrnd(var_params.mean, var_params.sd);
    
   % Draw assignments 
   MCsamples{l}.assignments = cell([N 1]);
   for ni = 1:N
    MCsamples{l}.assignments{ni} = postsynaptic_events{ni}.event_times;
    event_counts = size(postsynaptic_events{ni}.event_times,2);
    if event_counts > 0
        for event_i = 1:event_counts
            prob_tmp = [1; exp(var_params.phi{ni}(:,event_i))];
            prob_tmp = transpose(prob_tmp./sum(prob_tmp));
            assign_tmp =mnrnd(1,prob_tmp);
         MCsamples{l}.assignments{ni}(event_i) =  find(assign_tmp)-1;
        end  
    end
   end
end

%% Specify the prior distribution 
prior_parameters = struct('prob',0.2,'mu0',5,'sigma0',1);
%% Evaluate the stochastic gradient
Ngrid = 100;

fgrid = [1:Ngrid].*Tmax./Ngrid;
ftmp = zeros(Ngrid,1);
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
                        loglklh_p(cell_i) = loglklh_p(cell_i)+ log(ftmp(max(find(fgrid < event_times_tmp(event_i))))) + log(size_weights(event_i));
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
    
    % log q(a)
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
    
    
end

