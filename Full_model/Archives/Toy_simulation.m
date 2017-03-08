% Generate data from a simple toy example:
rng(12242,'twister');
%% Draw connectivity and event sizes 
% 
p=10;
conn_prob = 0.5;
gamma = unifrnd(0,1,p,1)<conn_prob;

size_sigma = 1;
size_mu = unifrnd(4,10,p,1); 

%% Draw evoked intensity for every cell in each trial 
N=500; % number of trials
induced_intensity = cell([N 1]);
K = 4;
for ni = 1:N
    induced_intensity{ni} = unifrnd(0.5,1,p,1); 
    % Make some zero
    induced_intensity{ni}(randsample(p,K)) = 0;
end

spontaneous_intensity = 0.2;
spontaneous_mu = 6;

%% Generate data 
% Note: the intensity is set to zero in the second after an event 
n0 = 100; %number of events drawn from dominating meaure.
Tmax=2;

% Initialize the storage 
presynaptic_events = cell([N p]);
postsynaptic_events = cell([N 1]);
for ni = 1:N
    dominating_measure = max(induced_intensity{ni});
    % Initialize the postsynaptic events with spontaneous events:
    postsynaptic_events{ni, 1} = struct('event_times',[],'event_size',[],'assignments',[]);
    
    event_candidate = cumsum(exprnd( 1/spontaneous_intensity, [n0 1]));
    postsynaptic_events{ni, 1}.event_times = transpose(event_candidate(event_candidate<Tmax));
    if sum(event_candidate<Tmax) >0 
        postsynaptic_events{ni, 1}.event_size = normrnd(spontaneous_mu,size_sigma,[1 sum(event_candidate<Tmax)]);
        postsynaptic_events{ni, 1}.assignments = zeros(1,sum(event_candidate<Tmax)); 
    end
    
    for pi = 1:p
        intensity_tmp = induced_intensity{ni}(pi);
        
        event_candidate = exprnd(1/dominating_measure, [n0 1]);
        
        event_tmp =[];
        size_tmp =[];
        
        event_counter = 0;
        t_now = 0;
        t_last = -1;
        while t_now < Tmax
            event_counter = event_counter+1;
            t_now = t_now + event_candidate(event_counter);
            if t_now < t_last + 1
                intensity_now = 0;
            else
                intensity_now= intensity_tmp;
            end
            
            if t_now < Tmax
                ac_prop = unifrnd(0,1);
                if ac_prop < intensity_now/dominating_measure
                    event_tmp = [event_tmp t_now];
                    w_now = normrnd(size_mu(pi),size_sigma);
                    size_tmp = [size_tmp w_now];
                    t_last = t_now;
                end
            end   
        end
        presynaptic_events{ni, pi} = struct('event_times',event_tmp,'event_size',size_tmp);
        if gamma(pi) == 1 
            postsynaptic_events{ni, 1}.event_times = [postsynaptic_events{ni, 1}.event_times event_tmp];
            postsynaptic_events{ni, 1}.event_size = [postsynaptic_events{ni, 1}.event_size size_tmp];
            if size(size_tmp,2)>0
            postsynaptic_events{ni, 1}.assignments = [postsynaptic_events{ni, 1}.assignments pi.*ones(size(size_tmp))];
            end
        end
        
    end
    
    % sort the events
    %postsynaptic_events{ni, 1}.event_times
	[events_this_trial events_order]= sort(postsynaptic_events{ni, 1}.event_times); 
    postsynaptic_events{ni, 1}.event_times = events_this_trial;
    postsynaptic_events{ni, 1}.event_size = postsynaptic_events{ni, 1}.event_size(events_order);
    postsynaptic_events{ni, 1}.assignments = postsynaptic_events{ni, 1}.assignments(events_order);
    
            
end

