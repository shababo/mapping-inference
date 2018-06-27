function [Kcor]=get_kernel_cor(X1,X2,tau)

          
    nsq1=sum(X1.^2,2);
    nsq2=sum(X2.^2,2);
    Kcor = nsq1*ones(1,length(nsq2)) +  ones(length(nsq1),1)*nsq2';
    Kcor=bsxfun(@minus,Kcor,(2*X1)*X2.');
    Kcor=exp(-Kcor/tau);
    