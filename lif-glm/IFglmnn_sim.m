%%
load('data\current_template.mat');
I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];

%% simulation

figure;suptitle('Simulation (by location)');
for locid=1:25
dt=1; %time step ms
t_end=75;t_StimStart=5;t_StimEnd=15;
V_th=-25; %spike threshold [mV]
E_L=-28; %resting membrane potential [mV]
V_reset=-85; %value to reset voltage to after a spike [mV]
g=0.1; %membrane time constant [ms]
k=0.025-0.0005*locid; %membrane resistance [MOhm]
k_vec(locid)=k;

stoc_mu=0;stoc_sigma=0.3;

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 
V_vect=zeros(1,length(t_vect));V_plot_vect=zeros(1,length(t_vect));
PlotNum=0;
I_Stim_vect=[100 100 100 100 100 50 50 50 50 50 25 25 25 25 25];
spTrain=zeros(t_end,length(I_Stim_vect));

for I_Stim=I_Stim_vect; %loop over different I_Stim values
    PlotNum=PlotNum+1;
    i=1; %index denoting which element of V is being assigned
    
    V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
    
    %%%%square current
%     I_e_vect=zeros(1,t_StimStart/dt); %portion of I_e_vect from t=0 to t_StimStart
%     I_e_vect=[I_e_vect I_Stim*ones(1,1+((t_StimEnd-t_StimStart)/dt))];
%     I_e_vect=[I_e_vect zeros(1,(t_end-t_StimEnd)/dt)];

    %%%%chi-sq shape current
    I_e_vect=[0;I_e(:,PlotNum)];
    
    I_e_vect_mat(:,PlotNum)=I_e_vect;
    
    NumSpikes=0; %holds number of spikes that have occurred
    for t=dt:dt:t_end %loop through values of t in steps of df ms        
%         V_inf = E_L + I_e_vect(i)*R_m; %value that V_vect is exponentially
%                                          %decaying towards at this time step
%         V_vect(i+1) = V_inf + (V_vect(i)-V_inf)*exp(-dt/tau);
%         V_vect(i+1) = V_vect(i) + (E_L-V_vect(i) + I_e_vect(i)*R_m)/tau*dt + sqrt(dt)*normrnd(stoc_mu,stoc_sigma);
        
        V_vect(i+1) = V_vect(i) + ((E_L-V_vect(i))*g + I_e_vect(i)*k)*dt + sqrt(dt)*normrnd(stoc_mu,stoc_sigma); 
        
        %%%% if statement below says what to do if voltage crosses threshold
        if (V_vect(i+1)>V_th) %cell spiked
            V_vect(i+1)=V_reset; %set voltage back to V_reset
            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spTrain(i,PlotNum)=1;
        end
 
        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
    end    
end

spTrain_byLoc{locid}=spTrain;
Ie_byLoc{locid}=I_e_vect_mat(1:end-1,:);

subplot(5,5,locid)
imagesc(spTrain');colormap(flipud(gray));box off;
set(gca,'YTick',5:5:15,'YTickLabel',5:5:15,'XTick',0:25:75,'XTickLabel',0:25:75,'FontSize',12);
hold on
plot([t_StimStart t_StimStart],[0 20],'k:');
plot([t_StimEnd t_StimEnd],[0 20],'k:');
plot([0 75],[5.5 5.5],'k-');
plot([0 75],[10.5 10.5],'k-');
plot([0 75],[15.5 15.5],'k-');
hold off
% xlabel('Time (ms)');
if locid==1
    ylabel('25mV 50mV 100mV');
end

[r100,n100]=max(spTrain(:,1:5)~=0,[],1);real_latency1st(1,locid)=mean(n100.*r100);
[r50,n50]=max(spTrain(:,6:10)~=0,[],1);real_latency1st(2,locid)=mean(n50.*r50);
[r25,n25]=max(spTrain(:,11:15)~=0,[],1);real_latency1st(3,locid)=mean(n25.*r25);
    
real_meanSp(1,locid)=round(mean(sum(spTrain(:,1:5),1))*100)/100;
real_stdSp(1,locid)=round(std(sum(spTrain(:,1:5),1))*100)/100;  
real_meanSp(2,locid)=round(mean(sum(spTrain(:,5:10),1))*100)/100;
real_stdSp(2,locid)=round(std(sum(spTrain(:,5:10),1))*100)/100;    
real_meanSp(3,locid)=round(mean(sum(spTrain(:,11:15),1))*100)/100;
real_stdSp(3,locid)=round(std(sum(spTrain(:,11:15),1))*100)/100;
end

simTrain_byStim=[];
for jj=1:15
    for ii=1:25
    simTrain_byStim{jj}(:,ii)=spTrain_byLoc{ii}(:,jj);
    end
end
figure;suptitle('Stimulation (by amplitude)');
for ii=1:15
    subplot(3,5,ii);
    imagesc(simTrain_byStim{ii}');box off;
    set(gca,'FontSize',12,'YTick',5:5:25,'YTickLabel',5:5:25);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    hold off
    if ii==1
        ylabel('Stim=100mV');
    elseif ii==6
        ylabel('Stim=50mV');
    elseif ii==11
        ylabel('Stim=25mV');
    end
    colormap(flipud(gray));
    xlabel('Time (ms)');
end

%% glm fit prep

%%%%%location dummy variable ind_allLoc
spTrain_allLoc=[];
Ie_allLoc=[];
ind1st_allLoc=[];
ind_allLoc=zeros(size(spTrain,1)*size(spTrain,2),25);
for ii=1:25
    spTrain=spTrain_byLoc{ii};
    ind1st=zeros(size(spTrain));
    for j=1:size(spTrain,2)
        sp1st=find(spTrain(:,j),1);
        ind1st(sp1st+1:end,j)=1; %%%%indicator for >1st spike
    end
    ind1st_allLoc{ii}=ind1st;
    spTrain_allLoc=[spTrain_allLoc;spTrain(:)];
    I_e=Ie_byLoc{ii};
    Ie_allLoc=[Ie_allLoc;I_e(:)];
    nii=size(spTrain,1)*size(spTrain,2);
    ind_allLoc((ii-1)*nii+1:ii*nii,ii)=1;
end

trainM=[];I_eg=[];ind1st=[];
for ii=1:25
    trainM=[trainM spTrain_byLoc{ii}];
    I_eg=[I_eg Ie_byLoc{ii}];
    ind1st=[ind1st ind1st_allLoc{ii}];
end


%% prepping design matrix

I_eg=repmat(I_e,1,25);
g=0.1;
[expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
nloc_sp=25;

expg_k_loc0=zeros(size(trainM(:),1),nloc_sp);
for ii=1:nloc_sp
    indloc=zeros(size(trainM,1)*size(trainM,2),1);
    indloc((ii-1)*nii+1:ii*nii)=1;
    expg_k_loc0(:,ii)=expg_k(:).*indloc;
end
% loc_min=find(locVec_sp==find(locVec_spC==min(locVec_spC(find(locVec_spC))),1));
loc_min=25;
expg_k_loc=expg_k_loc0(:,[1:loc_min-1 loc_min+1:end]);
design1=[ones(size(expg_k_loc,1),1) ind1st(:).*expg_Vreset(:)];
design2=[expg_k(:) expg_k_loc];
Y=trainM(:);

%% non-negative glm fit; call logL_LIFglmnn.m
design = [design1 design2];
num_params=size(design,2);
init_vals=zeros(1,num_params);
x0 = init_vals; x0(1)=-5;x0(2)=-50;
obj_fun = @(x) logL_LIFglmnn( x, design, Y );
options = optimset('Algorithm','interior-point','Display','iter','GradObj','on','Diagnostics','on','UseParallel','always'); %
                          
ub = Inf*ones(size(x0));                            
lb = -Inf*ones(size(x0));
lb(3:num_params) = -1e-6;
ub(3:num_params) = 1e6;

betahat_conv_allLoc = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);

%%
gain_vec=zeros(nloc_sp,1);
gain_vec(loc_min)=betahat_conv_allLoc(3);
loc_mmVec=[1:loc_min-1 loc_min+1:nloc_sp];
for i=1:nloc_sp-1
    gain_vec(loc_mmVec(i))=betahat_conv_allLoc(3)+betahat_conv_allLoc(3+i);
end

%% prediction
locVec=1:25;
figure;suptitle('Prediction (by amplitude)');
simTrainMat=zeros(75,9);
for jj=1:15
for ii=1:25
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    betahat_conv(3)=gain_vec(ii);

%% prediction
I_eg_fit=I_e(:,jj);
ntrial=1;
ntime=length(I_eg_fit);
simTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);step=zeros(ntime,1);
    j=1;lastSpT=0;numSp=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
    
    %%%%link function f(x) = 1+x, x>0; = exp(x), x<0.
    
    step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
    if step0>0
        step(j)=step0+1;
    else
        step(j)=exp(step0);
    end
    
    lambdaMat(j,tr)=step(j);
    
while (j<=ntime)
    if sum(step)>tao
        step(1:j)=0;
        simTrain(j,tr)=1; %spike
        tao=exprnd(1);
        lastSpT=j;numSp=numSp+1;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            
            step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
            if step0>0
                step(j)=step0+1;
            else
                step(j)=exp(step0);
            end
            
            lambdaMat(j,tr)=step(j);
        end
    else
        simTrain(j,tr)=0;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            
            if numSp==0
                step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            else
                step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            end
            
            lambdaMat(j,tr)=step(j);
        end
    end
end
end

simTrainMat(:,ii)=simTrain;
end
simTrainAll{jj}=simTrainMat;
subplot(3,5,jj);
imagesc(simTrainMat');
if jj==1
    ylabel('Stim=100mV');
elseif jj==6
    ylabel('Stim=50mV');
elseif jj==11
    ylabel('Stim=25mV');
end
set(gca,'YTick',5:5:25,'YTickLabel',5:5:25,'XTick',0:25:75,'XTickLabel',0:25:75,'FontSize',12);
hold on
plot([5 5],[0 122],'k:');
plot([15 15],[0 122],'k:');
hold off
colormap(flipud(gray));
xlabel('Time (ms)');
box off;
end

simTrain_byLoc=[];
for jj=1:15
    for ii=1:25
    simTrain_byLoc{ii}(:,jj)=simTrainAll{jj}(:,ii);
    end
end
figure;suptitle('Prediction (by location)');
for ii=1:25
    subplot(5,5,ii);
    imagesc(simTrain_byLoc{ii}');box off;
    set(gca,'YTick',5:5:15,'YTickLabel',5:5:15,'XTick',0:25:75,'XTickLabel',0:25:75,'FontSize',12);
    % title(['Loc=' num2str(locVec(ii))]);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    plot([0 75],[5.5 5.5],'k-');
    plot([0 75],[10.5 10.5],'k-');
    plot([0 75],[15.5 15.5],'k-');
    hold off
    if ii==1
        ylabel('25mV 50mV 100mV');
    end
    set(gca,'FontSize',12);
    colormap(flipud(gray));
    % xlabel('Time (ms)');
end

