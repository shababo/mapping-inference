function spikes = sim_glm(params, stims_x, stims_t, dt)

spikes = zeros(size(stims_t));
total_filter_out = params.baseline * ones(size(stims_t));
stim_filter_out = zeros(size(stims_t));

for i = 1:size(stims_t,1)
    
    stim_filter_conv_out = conv(stims_t(i,:),params.stim_filter);
    stim_filter_conv_out = stim_filter_conv_out(1:size(stims_t,2));
        
    total_filter_out(i,:) = total_filter_out(i,:) + ...
        params.spatial_footprint(stims_x(i,1)) * stim_filter_conv_out;
    
    stim_filter_out(i,:) = total_filter_out(i,:);
    
    % --------------- Set up simulation dynamics variables ---------------
    num_spikes = 0; % number of spikes
    current_t = 1; % current time bin
    next_spike_time_rescale = exprnd(1);  % time of next spike (in rescaled time)
    rprev = 0;  % Integrated rescaled time up to current point
    nbinsPerEval = 5;

    % --------------- Run dynamics ---------------------------------------
    while current_t <= size(stims_t,2)
        iinxt = current_t:min(current_t+nbinsPerEval-1,size(stims_t,2));
        rrnxt = exp(total_filter_out(i,iinxt))*dt; % Cond Intensity
        rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
        if (next_spike_time_rescale >= rrcum(end)) % No spike in this window
            current_t = iinxt(end)+1;
            rprev = rrcum(end);
        else   % Spike!
            ispk = iinxt(find(rrcum>=next_spike_time_rescale, 1, 'first')); % time bin where spike occurred
            num_spikes = num_spikes+1;
            spikes(i,ispk) = 1; % spike time
            mxi = min(size(stims_t,2), ispk+length(params.hist_filter)); % max time affected by post-spike kernel
            iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
            if ~isempty(iiPostSpk)
                total_filter_out(i,iiPostSpk) = ...
                    total_filter_out(i,iiPostSpk)+params.hist_filter(1:mxi-ispk);
            end
            next_spike_time_rescale = exprnd(1);  % draw next spike time
            rprev = 0; % reset integrated intensity
            current_t = ispk+1;  % Move to next bin
            % --  Update # of samples per iter ---
            muISI = current_t/num_spikes;
            nbinsPerEval = max(20, round(1.5*muISI)); 
        end
    end

end

assignin('base','total_filter_out',total_filter_out)
assignin('base','stim_filter_out',stim_filter_out)

