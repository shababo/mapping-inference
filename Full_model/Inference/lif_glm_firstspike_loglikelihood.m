function [loss] = lif_glm_firstspike_loglikelihood(x, responses, v_trace,background_rate,v_th_known,linkfunc,loss_type)

% linkfunc{1}: link
% linkfunc{2}: dlink
%%
% x: gain


% Calculate the firing probability
n_trial = size(v_trace,1);
n_grid = size(v_trace,2);
prob_firing = zeros(n_trial,n_grid);
for i_trial =1:n_trial
    for i_grid = 1:n_grid
        prob_firing(i_trial,i_grid)=linkfunc{3}(x(1)*v_trace(i_trial,i_grid)-v_th_known);
    end
end
   
if loss_type == 1
    
% 
prob_first_spike = zeros(n_trial,n_grid);

for i_trial =1:n_trial
    not_spike_prob=1;
    prob_first_spike(i_trial,1)=prob_firing(i_trial,1);
    for i_grid = 2:n_grid
    not_spike_prob = not_spike_prob*(1-prob_firing(i_trial,i_grid -1));
         prob_first_spike(i_trial, i_grid ) =not_spike_prob*prob_firing(i_trial,i_grid );
    end
end
   
least_square=zeros(n_trial,1);
for i_trial =1:n_trial
    least_square(i_trial)= ...
        sum( (background_rate+x(2)*prob_first_spike(i_trial,:)).^2);
    events=find(responses(i_trial,:)>0);
    if length(events)>0
        for i = 1:length(events)
            least_square(i_trial)= ...
                least_square(i_trial)-2*responses(i_trial,events(i))*(  x(2)*...
                prob_first_spike(i_trial,events(i))+background_rate);
        end
    end
end
loss= sum(least_square);

else

% the probability with background intensity:
loglklh=zeros(n_trial,1);
for i_trial =1:n_trial
    
    if max(0,(1-sum(responses(i_trial,:))))>0
    loglklh(i_trial)= ...
        (-background_rate*n_grid+log( 1-x(2)+x(2)*exp(sum( log(1- prob_firing(i_trial,:)) )) ))...
        *max(0,(1-sum(responses(i_trial,:))));
    end %otherwise it's zero
    events=find(responses(i_trial,:)>0);
    if length(events)>0
        for i = 1:length(events)
            if events(i)==1
                temp= log(1);
            else
                temp =sum( log(1- prob_firing(i_trial, 1:(events(i)-1) )));
            end
            loglklh(i_trial)= ...
                loglklh(i_trial)+responses(i_trial,events(i))*( log(  x(2)*exp(temp)*...
                prob_firing(i_trial,events(i))+background_rate));
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
