function [expg_Vreset,expg_EL,expg_k] = gconv(I_eg,trainM,g)

    ntime=size(I_eg,1);ntrial=size(I_eg,2);

    t_elapse=zeros(ntime,ntrial);
    expg_Vreset=zeros(ntime,ntrial);
    expg_EL=zeros(ntime,ntrial);
    expg_k=zeros(ntime,ntrial);

    for tr=1:ntrial
        if isempty(find(trainM(:,tr)))==1
            t_elapse=1:ntime;
            expg_Vreset(:,tr)=exp(-g.*t_elapse);
            expg_EL(:,tr)=1-exp(-g.*t_elapse);
            for t=1:ntime
                expg_k(t,tr)=exp(-g.*(t-[1:t]))*I_eg(1:t,tr);
            end
        elseif isempty(find(trainM(:,tr)))==0
            te0=find(trainM(:,tr));
            te1=[0;te0;ntime];
            for i=1:length(te0)+1
                t_elapse=1:(te1(i+1)-te1(i));
                expg_Vreset(te1(i)+1:te1(i+1),tr)=exp(-g.*t_elapse);
                expg_EL(te1(i)+1:te1(i+1),tr)=1-exp(-g.*t_elapse);
                for t=te1(i)+1:te1(i+1)
                    expg_k(t,tr)=exp(-g.*(t-[te1(i)+1:t]))*I_eg(te1(i)+1:t,tr);
                end
            end
        end
    end

end