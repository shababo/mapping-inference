function [ logL, grad ] = logL_LIFglmnn( beta, design, Y , link)
% Assuming that the spatial mark is know
% Therefore the third column of design matrix becomes the offset

num_params=size(beta,2);
num_T=size(design,1);


if link ==  1 %logit
    
    %%%% log-likelihood
    prob=zeros(num_T,1);
    logexb=zeros(num_T,1);logexb0=zeros(num_T,1);
    for i=1:num_T
        step=sum(beta.*design(i,1:num_params))+design(i,num_params+1);
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
    
    
elseif  link == 2 %log
    %
    %%%% log-likelihood
    step=zeros(num_T,1);
    lambda=zeros(num_T,1);
    for i=1:num_T
        step(i)=sum(beta.*design(i,1:num_params))+design(i,num_params+1);
        lambda(i)=exp(step(i));
    end
    logL=-sum(Y.*step-lambda);
    
    %%%% score function: gradient of logL
    grad=zeros(num_params,1);
    for j=1:num_params
        grad(j)=-sum(design(:,j).*(Y-lambda));
    end
    
end
%% link = linear rectifier

%%%% log-likelihood
% step=zeros(num_T,1);
% lambda=zeros(num_T,1);
% for i=1:num_T
%     step(i)=sum(beta.*design(i,1:3))+design(i,4);
%     lambda(i)=min(0.99,max(step(i),0.01));
% end
% logL= -sum(Y.*(log(lambda)) + (1-Y).*(log(1-lambda)));
%
% %%%% score function: gradient of logL
% grad=zeros(num_params,1);
% for j=1:num_params
%     grad(j)=-sum(design(:,j).* ((Y./lambda) + ((1-Y)./(1-lambda)) ));
%
% end


end
