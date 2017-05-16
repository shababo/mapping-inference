% independent copy of the gconv.m in ../lif-glm
function covs = gconv_v2(I_eg,stims_ind,trainM,g)

    ntime=size(I_eg,1);ntrial=size(I_eg,2);

    t_elapse=zeros(ntime,ntrial);
%     expg_Vreset=zeros(ntime,ntrial);
%     expg_EL=zeros(ntime,ntrial);
%     expg_k=zeros(ntime,ntrial);
    
    covs = zeros(ntime,2 + max(stims_ind),ntrial);
    covs(:,1,:) = -1;

    for tr=1:ntrial
        % no spikes
        if isempty(find(trainM(:,tr)))==1
            t_elapse=1:ntime;
            for t=1:ntime
                covs(t,2+stims_ind(tr),tr)=exp(-g.*(t-[1:t]))*I_eg(1:t,tr);
            end
        % spikes    
        elseif isempty(find(trainM(:,tr)))==0
            te0=find(trainM(:,tr));
            te1=[0;te0;ntime];
            
            for i=1:length(te0)+1
                if i ~= 1
                    t_elapse=1:(te1(i+1)-te1(i));
                    covs(te1(i)+1:te1(i+1),2,tr)=exp(-g.*t_elapse);
                end
                for t=te1(i)+1:te1(i+1)
                    covs(t,2+stims_ind(tr),tr)=exp(-g.*(t-[te1(i)+1:t]))*I_eg(te1(i)+1:t,tr);
                end
            end
        end
    end

end