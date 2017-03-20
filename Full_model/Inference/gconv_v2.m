function [cov_prior_first,cov_post_first,cov_current] = gconv_v2(I_eg,trainM,g,marks)

ntime=size(I_eg,1);ntrial=size(I_eg,2);
%expg_Vreset=zeros(ntime,ntrial);
%expg_EL=zeros(ntime,ntrial);
%expg_k=zeros(ntime,ntrial);
cov_prior_first=zeros(ntime,ntrial); % covariate before the first event
cov_post_first=zeros(ntime,ntrial); % covariate after the first event
cov_current=zeros(ntime,ntrial); % covariate for the current

for tr=1:ntrial
    if isempty(find(trainM(:,tr)))==1 % if there are no events
        t_elapse=1:ntime;
        cov_prior_first(:,tr)=exp(-g.*t_elapse);
        cov_post_first(:,tr) = 0;
        for t=1:ntime
            cov_current(t,tr)=exp(-g.*(t-[1:t]))*I_eg(1:t,tr);
            cov_prior_first(t,tr) = cov_prior_first(t,tr)+ g*sum(exp(-g.*(t-[1:t])));
        end
        cov_current(:,tr) = marks(tr).*cov_current(:,tr);
    elseif isempty(find(trainM(:,tr)))==0 % if there are more than one event
        te0=find(trainM(:,tr));
        te1=[0;te0;ntime];
        
        % Before the first event:
        t_elapse=1:min(te0);
        cov_prior_first(1:min(te0),tr)= exp(-g.*t_elapse);
        cov_post_first(1:min(te0),tr) = 0;
        for t=1:min(te0)
            cov_current(t,tr)= exp(-g.*(t-[1:t]))*I_eg(1:t,tr);
            cov_prior_first(t,tr) = cov_prior_first(t,tr)+ g*sum(exp(-g.*(t-[1:t])));
        end
        
        % After the first event
        for i= (1:length(te0))+1
            t_elapse=(te1(i)+1):te1(i+1);
            cov_prior_first(t_elapse,tr)= 0;
            cov_post_first(t_elapse,tr) = exp(-g.* (t_elapse-te1(i) ) );
            for t= (te1(i)+1):te1(i+1)
                cov_current(t,tr)=exp(-g.*(t-[te1(i)+1:t]))*I_eg(te1(i)+1:t,tr);
                cov_prior_first(t,tr) = cov_prior_first(t,tr)+ g*sum(exp(-g.*(t-[te1(i)+1:t])));
            end
            
        end
        cov_current(:,tr) = marks(tr).*cov_current(:,tr);
    end
end

end