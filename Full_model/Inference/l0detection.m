function [output] =l0detection(...
   y_trace,decay_rate,l0penalty,varargin)


% y_trace=mpp(i).trace;
% decay_rate = 0.98; %gamma in the paper
% l0penalty=16; %lambda in the paper 


if ~isempty(varargin) && ~isempty(varargin{1})
    
else
    
end

%% Apply the l0 spike detection:

F_zero= -l0penalty;
changepoints=cell(length(y_trace),1);
candidatesets=cell(length(y_trace),1);

F_value=zeros(length(y_trace),1);
F_candidates=zeros(length(y_trace),1);
D_value=zeros(length(y_trace),1);
C_value=zeros(length(y_trace),1);
changepoints{1}=1;
candidatesets{1}=0;
for s = 1:length(y_trace)
    F_candidates=zeros(length(candidatesets{s}),1);
   
    for tau_index = 1:length(candidatesets{s})
        tau=candidatesets{s}(tau_index)+1;
        decay_rate_vec=decay_rate.^(0:(s-tau));
        tmp_trace=y_trace(tau:s);
        C_value(tau)= sum(tmp_trace.*decay_rate_vec)/sum(decay_rate_vec.^2);
        D_value(tau)= sum(tmp_trace.^2)/2-C_value(tau)*sum(tmp_trace.*decay_rate_vec)+...
            C_value(tau)^2*sum(decay_rate_vec.^2)/2;
        if tau == 1
            F_candidates(tau_index)=F_zero+D_value(tau)+l0penalty;
        else
            F_candidates(tau_index)=F_value(tau-1)+D_value(tau)+l0penalty;
        end
    end
    [F_value(s), s_prime_index]=min(F_candidates);
    s_prime=candidatesets{s}(s_prime_index);
    candidatesets{s+1}=[candidatesets{s}((F_candidates(1:length(candidatesets{s})) -l0penalty)< F_value(s)); s];
    if s_prime ~= 0
    changepoints{s+1}=[changepoints{s_prime} s_prime];
    end
end
%%
%Calculate the estimated calcium concentrations
y_fit=zeros(length(y_trace),1);
for s = 1:length(y_trace)
    prior_changepoints=find(changepoints{end}< s);
    if ~isempty(prior_changepoints)
       last_changepoint=changepoints{end}(max(prior_changepoints));
       if max(prior_changepoints)==length(changepoints{end})
           next_changepoint = length(y_trace);
       else
           next_changepoint=changepoints{end}(max(prior_changepoints)+1);
       end
       if s~=(last_changepoint+1)
            y_fit(s)=y_fit(s-1)*decay_rate;
       else
           decay_rate_vec=decay_rate.^(0:(next_changepoint-last_changepoint-1));
           tmp_trace=y_trace((last_changepoint+1):next_changepoint);
           y_fit(s)= sum(tmp_trace.*decay_rate_vec)/sum(decay_rate_vec.^2);
       end
    else
        
    end
end
output=struct;
output.trace=y_trace;
output.fit=y_fit;
if ~isempty(changepoints{end})
output.changepoints=changepoints{end}+1;
else
    output.changepoints=[];
end

