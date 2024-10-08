function [ logL, grad ] = logL_LIFglmnn( beta, design, Y, linkID)


num_params=size(design,2);
num_T=size(design,1);


if linkID==1
%% link = identity: for fitting voltage
%%%% log-likelihood
lambda=zeros(num_T,1);
for i=1:num_T
    lambda(i)=sum(beta.*design(i,:));
end
logL=sum((Y-lambda).^2);

%%%% score function: gradient of logL
grad=zeros(num_params,1);
for j=1:num_params
    grad(j)=-sum(design(:,j).*(Y-lambda));
end


elseif linkID==2
%% link = log: hard spike assignment
%%%% log-likelihood
step=zeros(num_T,1);
lambda=zeros(num_T,1);
for i=1:num_T
    step(i)=sum(beta.*design(i,:));
    lambda(i)=exp(step(i));
end
logL=-sum(Y.*step-lambda);

%%%% score function: gradient of logL
grad=zeros(num_params,1);
for j=1:num_params
    grad(j)=-sum(design(:,j).*(Y-lambda));
end
    
% elseif linkID==3
% %% link = soft rectifier: hard spike assignment


elseif linkID==4
%% link = logit: soft spike assignment
%%%% log-likelihood
prob=zeros(num_T,1);
logexb=zeros(num_T,1);logexb0=zeros(num_T,1);
for i=1:num_T
    step=sum(beta.*design(i,:));
    prob(i)=exp(step)/(1+exp(step));
    logexb(i)=log(prob(i));
    logexb0(i)=log(1-prob(i));
end
logL= -sum(Y.*logexb + (1-Y).*logexb0);

%%%% score function: gradient of logL
grad=zeros(num_params,1);
for j=1:num_params
    grad(j)=-sum(design(:,j).*(Y-prob));
end

end

end
