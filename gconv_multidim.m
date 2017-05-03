function [expg_Vreset,expg_EL,expg_k] = gconv_multidim(input_current,spikes,g)

    ntime=size(input_current,1);ntrial=size(input_current,2);

    t_elapse=zeros(ntime,ntrial);
    expg_Vreset=zeros(ntime,ntrial);
    expg_EL=zeros(ntime,ntrial);
    expg_k=zeros(ntime,ntrial);
    expg_dk=zeros(ntime,ntrial);

    for tr=1:ntrial
        if isempty(find(spikes(:,tr)))==1
            t_elapse=1:ntime;
%             expg_Vreset(:,tr)=zeros(size(t_elapse));
            expg_EL(:,tr)=1-exp(-g.*t_elapse);
            for t=1:ntime
                expg_k(t,tr)=exp(-g.*(t-[1:t]))*input_current(1:t,tr);
%                 expg_dk(t,tr)=exp(-g.*(t-[1:t]))*d_input_current(1:t,tr);
            end
        elseif isempty(find(spikes(:,tr)))==0
            te0=find(spikes(:,tr));
            te1=[0;te0;ntime];
            for i=1:length(te0)+1
                t_elapse=1:(te1(i+1)-te1(i));
                if i ~= 1
                    expg_Vreset(te1(i)+1:te1(i+1),tr)=exp(-g.*t_elapse);
                    expg_EL(te1(i)+1:te1(i+1),tr)=1-exp(-g.*t_elapse);
                end
                for t=te1(i)+1:te1(i+1)
                    expg_k(t,tr)=exp(-g.*(t-[te1(i)+1:t]))*input_current(te1(i)+1:t,tr);
%                     expg_dk(t,tr)=exp(-g.*(t-[te1(i)+1:t]))*d_input_current(te1(i)+1:t,tr);
                end
            end
        end
    end

end

