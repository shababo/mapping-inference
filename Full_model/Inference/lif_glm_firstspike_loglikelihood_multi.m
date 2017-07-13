function [loss] = lif_glm_firstspike_loglikelihood_multi(x, responses, ...
    v_trace,background_rate,v_th_known,linkfunc,loss_type,delay_params)

% linkfunc{1}: link
% linkfunc{2}: dlink

% x: gain & gamma
n_cell=size(v_trace,1);
n_trial = size(v_trace,2);
n_grid = size(v_trace,3);

prob_firing=zeros(n_cell,n_trial,n_grid);
for i_cell = 1:n_cell
    for i_trial =1:n_trial
        for i_grid = 1:n_grid
            prob_firing(i_cell,i_trial,i_grid)=linkfunc{3}(x(2*(i_cell-1) +1)*v_trace(i_cell,i_trial,i_grid)-v_th_known(i_cell));
        end
    end
end

prob_first_spike=zeros(n_cell,n_trial,n_grid);
for i_cell = 1:n_cell
    for i_trial =1:n_trial
        not_spike_prob=1;
        prob_first_spike(i_cell,i_trial,1)=prob_firing(i_cell,i_trial,1);
        for i_grid = 2:n_grid
            not_spike_prob = not_spike_prob*(1-prob_firing(i_cell,i_trial,i_grid -1));
            prob_first_spike(i_cell,i_trial, i_grid) =not_spike_prob*prob_firing(i_cell,i_trial,i_grid );
        end
    end
end
% Convoluted with delays:
% Calculate the delay density on a grid

if delay_params.delayed == 0
    prob_first_spike_delayed = prob_first_spike;% do nothing 
else
    delay_prob = zeros(delay_params.n_grid+1,1);
    if delay_params.type == 1 % normal
        delay_prob = normpdf((0:delay_params.n_grid),delay_params.mean,...
            delay_params.std);
    elseif delay_params.type == 2 % gamma
        shape=(delay_params.mean^2)/(delay_params.std^2);
        %scale 
        scale = delay_params.mean/shape;
        delay_prob = gampdf((0:delay_params.n_grid),shape,scale);
    end
    % we approximate the probability with densities
    delay_prob = delay_prob/sum(delay_prob);
    min_delay = 0;
    max_delay = delay_params.n_grid;
    prob_first_spike_delayed = prob_first_spike;
    
    for i_cell = 1:n_cell
        for i_trial = 1:n_trial
            for i_grid = 1:n_grid
                
                idx_time = max(i_grid-max_delay,1): min(i_grid-min_delay,delay_params.n_grid);
                idx_delay = -( (min(idx_time)-i_grid) : (max(idx_time)-i_grid))+1;
                %M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
                temp=0;
                for i = 1:length(idx_time)
                    temp=temp+...
                        prob_first_spike(i_cell,i_trial,idx_time(i))*delay_prob(idx_delay(i));
                end
                prob_first_spike_delayed(i_cell,i_trial, i_grid) =temp;
                
            end
        end
    end
end

prob_no_spike = zeros(n_cell,n_trial);
for i_cell = 1:n_cell
    for i_trial =1:n_trial
        prob_no_spike(i_cell,i_trial)= 1-sum(prob_first_spike_delayed(i_cell,i_trial, :));
    end
end

%plot(reshape(prob_firing(i_cell,i_trial,:), [n_grid,1]))
if loss_type == 1
    least_square=zeros(n_trial,1);
    for i_trial =1:n_trial
        temp=zeros(1,1,n_grid);
        for i_cell = 1:n_cell
            temp=temp+x((i_cell-1)*2+2)*prob_first_spike(i_cell,i_trial,:);
        end
        least_square(i_trial)= ...
            sum( (temp+background_rate).^2);
        events=find(responses(i_trial,:)>0);
        if length(events)>0
            for i = 1:length(events)
                temp=background_rate;
                for i_cell = 1:n_cell
                    temp=temp+x((i_cell-1)*2+2)*prob_first_spike(i_cell,i_trial,events(i));
                end
                least_square(i_trial)= least_square(i_trial)-2*responses(i_trial,events(i))*temp;
            end
        end
    end
    loss= sum(least_square);
    
else
    
    % the probability with background intensity:
    loglklh=zeros(n_trial,1);
    for i_trial =1:n_trial
        if max(0,(1-sum(responses(i_trial,:))))>0
            temp=0;
            for i_cell= 1:n_cell
                %temp_prob=1-x((i_cell-1)*2+2)+x((i_cell-1)*2+2)*...
                %   exp(sum( log(1- prob_firing(i_cell,i_trial,:)) )) ;
                temp_prob=...
                    1-x((i_cell-1)*2+2)+x((i_cell-1)*2+2)*prob_no_spike(i_cell,i_trial);
                if temp_prob <1e-8
                    temp_prob=1e-8; % avoid singularity
                end
                temp=temp+log(temp_prob);
            end
            loglklh(i_trial)= (-background_rate*n_grid+temp)...
                *max(0,(1-sum(responses(i_trial,:))));
        end %otherwise it's zero
        
        events=find(responses(i_trial,:)>0);
        if length(events)>0
            for i = 1:length(events)
                temp2=background_rate;
                for i_cell = 1:n_cell
                    %                     if events(i)==1
                    %                         temp= log(1);
                    %                     else
                    %                         temp=sum( log(1- prob_firing(i_cell,i_trial, 1:(events(i)-1) )));
                    %                     end
                    %                     temp2=temp2+x((i_cell-1)*2+2)*exp(temp)*...
                    %                         prob_firing(i_cell,i_trial,events(i));
                    temp2=temp2+x((i_cell-1)*2+2)*prob_first_spike_delayed(i_cell,i_trial,events(i));
                    %                     prob_first_spike(i_cell,i_trial,events(i));
                end
                loglklh(i_trial)= ...
                    loglklh(i_trial)+responses(i_trial,events(i))*log(temp2);
            end
        end
    end
    loss= -sum(loglklh);
end


% Calculate the (expected) log-likelihood
% loglklh=zeros(n_trial,1);
% for i_trial =1:n_trial
%     loglklh(i_trial)= sum(log(1- prob_firing(i_trial,:))  )*max(0,(1-sum(responses(i_trial,:))));
%     events=find(responses(i_trial,:)>0);
%     if length(events)>0
%         for i = 1:length(events)
%             if events(i)==1
%                 temp= log(1);
%             else
%                 temp =sum( log(1- prob_firing(i_trial, 1:(events(i)-1) )));
%             end
%             loglklh(i_trial)= ...
%                 loglklh(i_trial)+responses(i_trial,events(i))*( temp+...
%                 log(prob_firing(i_trial,events(i))));
%         end
%     end
% end
% logL= -sum(loglklh);


% For gradient:

%
% % Calculate the partial derivtive w.r.t. x
% partial_x_prob = zeros(n_trial,n_grid);
% for i_trial =1:n_trial
%     for i_grid = 1:n_grid
%         partial_x_prob(i_trial,i_grid)=linkfunc{4}(x*v_trace(i_trial,i_grid)-v_th_known)*v_trace(i_trial,i_grid);
%     end
% end
%
% grad_samples=zeros(n_trial,1);
% for i_trial =1:n_trial
%     grad_samples(i_trial)=  sum(-partial_x_prob(i_trial,:)./(1- prob_firing(i_trial,:)))*max(0,(1-sum(responses(i_trial,:))));
%     events=find(responses(i_trial,:)>0);
%     if length(events)>0
%         for i = 1:length(events)
%             if events(i)==1
%                 temp = 0;
%             else
%                 temp=sum(-partial_x_prob(i_trial,1:(events(i)-1))./(1- prob_firing(i_trial,1:(events(i)-1))));
%             end
%             grad_samples(i_trial)= ...
%                 grad_samples(i_trial)+responses(i_trial,events(i))*...
%                 (temp ...
%                 +partial_x_prob(i_trial,events(i))/prob_firing(i_trial,events(i)));
%         end
%     end
% end
% grad= -sum(grad_samples);
%


% Calculate the firing probability
% n_trial = size(v_trace,1);
% n_grid = size(v_trace,2);
% prob_firing = zeros(n_trial,n_grid);
% max_time = zeros(n_trial,1);
% for i_trial =1:n_trial
%     for i_grid = 1:n_grid
%         prob_firing(i_trial,i_grid)=linkfunc{3}(x*v_trace(i_trial,i_grid)-v_th_known);
%         if sum(prob_firing(i_trial,1:i_grid)) < 1
%         max_time(i_trial) = i_grid;
%         end
%     end
% end
% % Calculate the (expected) log-likelihood
% loglklh=zeros(n_trial,1);
% for i_trial =1:n_trial
%     if (max_time(i_trial)+1> n_grid)
%         loglklh(i_trial)= (log(1- sum(prob_firing(i_trial,:)) )*max(0,(1-sum(responses(i_trial,1)))));
%     end
%     events=find(responses(i_trial,:)>0);
%     if length(events)>0
%         for i = 1:length(events)
%             if events(i) < max_time(i_trial)+1
%                 loglklh(i_trial)= ...
%                     loglklh(i_trial)+responses(i_trial,events(i))*( log(1- sum(prob_firing(i_trial,events(i)-1))) +...
%                     log(prob_firing(i_trial,events(i))));
%             end
%         end
%     end
% end
% logL= -sum(loglklh);

% For gradient:


% % Calculate the partial derivtive w.r.t. x
% partial_x_prob = zeros(n_trial,n_grid);
% for i_trial =1:n_trial
%     for i_grid = 1:n_grid
%         partial_x_prob(i_trial,i_grid)=linkfunc{4}(x*v_trace(i_trial,i_grid)-v_th_known)*v_trace(i_trial,i_grid);
%     end
% end
%
% grad_samples=zeros(n_trial,1);
% for i_trial =1:n_trial
%     if (max_time(i_trial)+1> n_grid)
%         grad_samples(i_trial)=  (-sum(partial_x_prob(i_trial,:)))/(1- sum(prob_firing(i_trial,:)) )*max(0,(1-sum(responses(i_trial,1))));
%     end
%     events=find(responses(i_trial,:)>0);
%     if length(events)>0
%         for i = 1:length(events)
%             if events(i) < max_time(i_trial)+1
%                 grad_samples(i_trial)= ...
%                     grad_samples(i_trial)+responses(i_trial,events(i))*...
%                     ((-sum(partial_x_prob(i_trial,1: (events(i)-1)) ))/(1- sum(prob_firing(i_trial, 1:(events(i)-1) )) )+...
%                     +partial_x_prob(i_trial,events(i))/prob_firing(i_trial,events(i)));
%             end
%         end
%     end
% end
% grad= -sum(grad_samples);
