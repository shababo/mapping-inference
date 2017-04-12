%%%% demo LIFglmnn.m on SOM data
%%%% call gconv.m and logL_LIFglmnn.m

load('som_lif_glm_data.mat');
load('som_lif_glm_voltage_data.mat');

%% choose cell 1
ind_e_all=[3 5]; 
ind_sp_all=[1 4];
ind_v_all=[1 3]; %intracellular voltage

%%
ind_sp=ind_sp_all(1);ind_e=ind_e_all(1);ind_v=ind_v_all(1);
locVec=1:25;nloc=length(locVec);

spTrain=[];spTrain_byLoc=[];
Ie_byLoc=[];Vm_byLoc=[];
for loc=locVec
    spTrain=all_data{ind_sp}{loc}';
    Ie=-repmat(all_data{ind_e}{7}(:,1:20:2000)'/min(all_data{ind_e}{7}(:,1:20:2000)),1,9)*-1; % do
    Vm=voltage_data{ind_v}{loc}(:,1:20:2000)';

    spTrain_byLoc{loc}=spTrain;
    Ie_byLoc{loc}=Ie;
    Vm_byLoc{loc}=Vm;  
end


%%
% linkID=1; %%% linear regression
% [ gain_vec_lm ] = LIFglmnn( spTrain_byLoc, Ie_byLoc, Vm_byLoc, linkID);

linkID=2; %%% glm
[ gain_vec_glm ] = LIFglmnn( spTrain_byLoc, Ie_byLoc, Vm_byLoc, linkID);



figure;plot(gain_vec_lm,'b*');hold on;
plot(gain_vec_glm,'rO');
hold off;
set(gca,'FontSize',16);
xlabel('Locations');ylabel('Spatial weights');

