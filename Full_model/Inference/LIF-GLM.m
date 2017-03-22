%% Fit the LIF-GLM model with soft-assignments 
% Updates: allow for more than two spikes, and allow for offsef in the GLM
% fit

%%%%%location dummy variable ind_allLoc

spTrain_allLoc=[];
Ie_allLoc=[];
ind1st_allLoc=[];
ind_allLoc=zeros(size(spTrain,1)*size(spTrain,2),25);
%%
% Modify the indicator function to allow for more than two spikes..
ii=5;

spTrain=spTrain_byLoc{ii};
ind1st=zeros(size(spTrain));
for j=1:size(spTrain,2)
    sp1st=find(spTrain(:,j),1);
    ind1st(sp1st+1:end,j)=1; %%%%indicator for >1st spike
end

ind1st_allLoc{ii}=ind1st;
spTrain_allLoc=[spTrain(:)];
I_e=Ie_byLoc{ii};
Ie_allLoc=[I_e(:)];

nii=size(spTrain,1)*size(spTrain,2); % Total number of observed time points?

ind_allLoc((ii-1)*nii+1:ii*nii,ii)=1; % some indicator function


trainM=[];I_eg=[];ind1st=[];

trainM=[spTrain_byLoc{ii}];
I_eg=[Ie_byLoc{ii}];
ind1st=[ind1st_allLoc{ii}];

%%
I_eg=repmat(I_e,1,1); % Stimuli size
marks = k_vec(ii).*ones(size(spTrain,2),1);
g=0.1;
[cov1,cov2,cov3]=gconv_v2(I_eg,trainM,g,marks); %temporally convolve paramters with g upto spike time

design=[ -ones(size(cov1(:))) cov1(:) cov2(:) cov3(:)];
%design=[ones(size(cov1(:))) cov2(:) cov3(:)];
Y=trainM(:);

%% non-negative glm fit; call logL_LIFglmnn.m
num_params=size(design,2)-1;
init_vals=zeros(1,num_params);

x0 = init_vals; x0(1)=-5;x0(2)=-20;x0(3)=-50;
obj_fun = @(x) logL_LIFglmnn_offset( x, design, Y );
options = optimset('Algorithm','interior-point','Display','iter','GradObj','on','Diagnostics','on','UseParallel','always'); %
                          
ub = Inf*ones(size(x0));                            
lb = -Inf*ones(size(x0));
lb(3:num_params) = -1e-6;
ub(3:num_params) = 1e6;

betahat_conv_allLoc = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);

%%
betahat_conv_allLoc-85