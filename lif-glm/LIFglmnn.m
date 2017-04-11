function [ gain_vec ] = LIFglmnn( sptrain_byLoc, Ie_byLoc, Vm_byLoc, linkID )
%%%%Ben:
% basically a function that takes in 
% spike times, stim info (power, location, current template), and some params.
% until we have a true shape model the locations are independent:
% a function that has no location info (and you input just one location's data),
% then you can just write a function around this that calls it for each location in the data
%%%%
%%% use the volatge (i.e. i-clamp) data to see how well our model is doing,
%%% but in the end it's the juxta we should fit
%%% the more "accurate" data is the juxta data

%% Input
%1. sptrain_byLoc: a cell array,
%   where each cell i-th is a TxN matrix of spike train at loc=i
%   and T=length of a trial and N=number of trial.
%2. Ie_byLoc: a cell array; corresponding v-clamp current.
%3. Vm_byLoc: a cell array; corresponding i-clamp voltage.
%4. linkID: choice of link function; 1=identity; 2=log; (3=soft rectifire;) 4=logit.

%% Output
% gain_vec: the vector of spatial weights (length = # of location)

%%
nloc=size(sptrain_byLoc,2);
loc_vec=1:nloc;

trainM=[];IeM=[];VmM=[];
for iloc=loc_vec
    trainM=[trainM sptrain_byLoc{iloc}];
    IeM=[IeM Ie_byLoc{iloc}]; % current known at every location
    VmM=[VmM Vm_byLoc{iloc}];
end

ind1st_allLoc=[];ind1st=[];spC=zeros(nloc,1);
for iloc=loc_vec
    sptrain=sptrain_byLoc{iloc};
    ind1st0=zeros(size(sptrain));
    for j=1:size(sptrain,2)
        sp1st=find(sptrain(:,j),1);ind1st0(sp1st+1:end,j)=1; %%%%indicator for >1st spike
    end
    spC(iloc)=sum(sptrain(:)); %%% total spike counts at each location
    ind1st_allLoc{iloc}=ind1st0;ind1st=[ind1st ind1st0];
end

%%%% truncate v-clamp voltage uptill first spike time
nii=size(sptrain,1)*size(sptrain,2);
Ie_loc0=zeros(size(trainM(:),1),nloc);
Ie_loc0_trunc=Ie_loc0;
for ii=1:nloc
    indloc=zeros(size(trainM,1)*size(trainM,2),1);
    indloc((ii-1)*nii+1:ii*nii)=1;
    Ie_loc0(:,ii)=IeM(:).*indloc;
    Ie_loc0_trunc(:,ii)=IeM(:).*indloc.*(1-ind1st(:));
end
loc_min=find(spC==0,1);
Ie_loc=Ie_loc0(:,[1:loc_min-1 loc_min+1:end]);
Ie_loc_trunc=Ie_loc0_trunc(:,[1:loc_min-1 loc_min+1:end]);



if linkID==1 %% lm non-negative
    IeM_trunc=IeM.*(1-ind1st);VmM_trunc=VmM.*(1-ind1st);
    design=[ones(size(IeM_trunc(:))) IeM_trunc(:) Ie_loc_trunc];Y=VmM_trunc(:);
    
    num_params=size(design,2);
    init_vals=zeros(1,num_params);
    x0 = init_vals;
    lb=-Inf*ones(size(x0));
    ub=Inf*ones(size(x0));                            
    lb(2:num_params)=-1e-6;
    ub(2:num_params)=1e6;
    
elseif (linkID==2 | linkID==3 |linkID==4) %% glm non-negative, convolve
    g=0.1;
    nii=size(sptrain,1)*size(sptrain,2);
    [expg_Vreset,expg_EL,expg_k]=gconv(IeM,trainM,g); %temporally convolve paramters with g upto spike time
    expg_k_loc0=zeros(size(trainM(:),1),nloc);
    for ii=1:nloc
        indloc=zeros(size(trainM,1)*size(trainM,2),1);
        indloc((ii-1)*nii+1:ii*nii)=1;
        expg_k_loc0(:,ii)=expg_k(:).*indloc;
    end
    expg_k_loc=expg_k_loc0(:,[1:loc_min-1 loc_min+1:end]);

    design1=[ones(size(expg_k_loc,1),1) ind1st(:).*expg_Vreset(:)];
    design2=[expg_k(:) expg_k_loc];
    design = [design1 design2];
    Y=trainM(:);
    
    num_params=size(design,2);
    init_vals=zeros(1,num_params);
    x0 = init_vals; x0(1)=-5;x0(2)=-50;
    ub = Inf*ones(size(x0));                            
    lb = -Inf*ones(size(x0));
    lb(3:num_params) = -1e-6;
    ub(3:num_params) = 1e6;
end

obj_fun = @(x) logL_LIFglmnn( x, design, Y, linkID);
options = optimset('Algorithm','interior-point','Display','iter','GradObj','on','Diagnostics','on','UseParallel','always'); %

betahat_conv_allLoc = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);

if linkID==1 %% lm non-negative
    gain_vec=zeros(nloc,1);
    gain_vec(loc_min)=betahat_conv_allLoc(2);
    loc_mmVec=[1:loc_min-1 loc_min+1:nloc];
    for i=1:nloc-1
        gain_vec(loc_mmVec(i))=betahat_conv_allLoc(2)+betahat_conv_allLoc(2+i);
    end

elseif (linkID==2 | linkID==3 |linkID==4) %% glm non-negative, convolve
    gain_vec=zeros(nloc,1);
    gain_vec(loc_min)=betahat_conv_allLoc(3);
    loc_mmVec=[1:loc_min-1 loc_min+1:nloc];
    for i=1:nloc-1
        gain_vec(loc_mmVec(i))=betahat_conv_allLoc(3)+betahat_conv_allLoc(3+i);
    end    
end


% figure; plot(gain_vec,'b*');set(gca,'FontSize',16);
% xlabel('Locations');ylabel('Spatial weights');


end

