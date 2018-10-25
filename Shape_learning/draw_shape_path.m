function [ ]=draw_shape_path(one_neuron, params,prior_params, prior_info)
% params.plot.shape_grid controls 
%% Record the path of fits
% % This is basically a reformat of the parameter_history
% parameter_path =struct;
% n_iterations=size(parameter_history,1);
% n_neurons=size(parameter_history,2);
% fldnames = fieldnames(parameter_history(1,1));
% 
% for i=1:length(fldnames)
%     i_field = fldnames{i};
%     if strcmp(i_field,'shapes')
%         parameter_path.(i_field)=struct;
%         parameter_path.(i_field).neurons(n_neurons)=struct;
%         for i_neuron = 1:n_neurons
%             tmp_mat = zeros(n_iterations,length(parameter_history(1,i_neuron).(i_field).mean));
%             tmp_mat_sig= tmp_mat;
%             for i_iteration=1:n_iterations
%                 tmp_mat(i_iteration,:) = parameter_history(i_iteration,i_neuron).(i_field).mean;
%                 tmp_mat_sig(i_iteration,:) = exp(parameter_history(i_iteration,i_neuron).(i_field).log_sigma);
%             end
%             parameter_path.(i_field).mean=tmp_mat; parameter_path.(i_field).sigma=tmp_mat_sig
%             parameter_path.(i_field).neurons(i_neuron).location=parameter_history(1,i_neuron).shapes.locations;
%             parameter_path.(i_field).neurons(i_neuron).bounds=parameter_history(1,i_neuron).shapes.bounds;
%         end
%     end
% end
%% Parameters:
% sigma_epsilon =  0.1;
x_grid_sample =-30:2:30;
y_grid_sample = -30:2:30;
  S=200;
n_grid =51;
this_z =0;
%% Record the true values
draw_truth=false;
if isfield(params,'truth')
    draw_truth=true;
end

%% Draw plots
%% True shape:
  % Calculate the shape values at the mesh grid (similar to visualize_3D_GP

this_z_indices=find(params.mesh_grid(:,3)==this_z);
true_shape_values=one_neuron.truth.shape(this_z_indices);
xy_locations=params.mesh_grid(this_z_indices,1:2);
x_grid = (1:n_grid)*(range(xy_locations(:,1))/n_grid) +min(xy_locations(:,1));
y_grid = (1:n_grid)*(range(xy_locations(:,2))/n_grid) +min(xy_locations(:,2));
% y_grid = (1:n_grid)*(80/n_grid) -40;
 cscale=   [min(true_shape_values) max(true_shape_values)];
[Xq,Yq] = meshgrid(x_grid,y_grid);
%             title(['z plane:' num2str(unique_z(i_z)) ])
%% Pick out the selected parameter sets
clear('param_list')
param_list(1)=prior_params;
for i =2:length(one_neuron.params)
    param_list(i)=one_neuron.params(i);
end

%% Draw samples from the prior shape
[X_tmp,Y_tmp] = meshgrid(x_grid_sample,y_grid_sample);
x_vec=reshape(X_tmp,[],1);y_vec=reshape(Y_tmp,[],1);
new_locations= zeros(length(x_vec),3);
new_locations(:,1)=x_vec;new_locations(:,2)=y_vec;
%% Draw samples from each of the chosen iterations:
posterior_samples=cell([S length(param_list)]);
post_mean_shape=cell([length(param_list) 1]);
clear('new_shape_params')
post_var_shape=cell([length(param_list) 1]);
for i=1:length(param_list)
    current_params=param_list(i);
    if params.plot.prediction_joint
        [new_shape_params(i)]=get_new_shape_joint(current_params,new_locations,prior_info);
%         for s=1:S
%             posterior_samples{s,i}=mvnrnd(new_shape_params(i).post_mean,new_shape_params(i).Sigma_tilde)';
%         end
    else
        [new_shape_params(i)]=get_new_shape_conditional(current_params,new_locations,prior_info);
%         for s=1:S
%             posterior_sample = draw_samples_from_var_dist(current_params);
%             post_new_sample=draw_samples_from_shape_conditional(new_shape_params(i),posterior_sample);
%             posterior_samples{s,i}=[posterior_sample.shapes; post_new_sample.shapes];
%           
%         end
    end
%     post_mean_shape{i}= [posterior_samples{1,i}]/S;
%     for s=2:S
%         post_mean_shape{i}= post_mean_shape{i}+[posterior_samples{s,i}]/S;
%     end
%     % calculate the posterior variance from the samples:
%     post_var_shape{i}=zeros(length(post_mean_shape{i}),1);
%         for i_loc = 1:length(post_mean_shape{i})
%        tmp=zeros(S,1);
%         for s=1:S
%         tmp(s)=   posterior_samples{s,i}(i_loc);
%         end
%         post_var_shape{i}(i_loc)=var(tmp);
%     end
tmp=current_params.shapes.mean;
tmp_mean = exp(tmp)./(1+exp(tmp));
 tmp_new=new_shape_params(i).mean_new + new_shape_params(i).Sigma_mean*(tmp_mean- new_shape_params(i).mean_old); %
post_mean_shape{i} =[tmp_mean;tmp_new];

% Conditional variance equation:
new_var=new_shape_params(i).Sigma_cond+ new_shape_params(i).Sigma_mean*current_params.shapes.Sigma_tilde*new_shape_params(i).Sigma_mean';
   post_var_shape{i} =[diag(current_params.shapes.Sigma_tilde);diag(new_var)];
end
%% Draw plots: 

if params.plot.do
 % Prior to truth:
 Vq = griddata(xy_locations(:,1),xy_locations(:,2),true_shape_values,Xq,Yq);

 cscale=   [min(true_shape_values) max(true_shape_values)];
 color_list= flip(copper(length(param_list)+1));
 figure(1)
for i=1:length(param_list)
 subplot(1,length(param_list)+1,i) 
Vq_fits = griddata(new_shape_params(i).locations(:,1),new_shape_params(i).locations(:,2),post_mean_shape{i},Xq,Yq);
imagesc(y_grid,x_grid,Vq_fits',cscale)
hold on;
%     caxis(cscale)

% hold on;
if i==1
    ylabel('x');xlabel('y');title(['Prior shape'])  
else
    scatter(new_shape_params(i).old_locations(:,2),new_shape_params(i).old_locations(:,1) ,'MarkerEdgeColor',color_list(i,:)) % stimulated locations

 ylabel('x');xlabel('y');title(['Batch ' num2str(i-1)])  
end
colorbar('eastoutside')
end

 subplot(1,length(param_list)+1,i+1) 
imagesc(y_grid,x_grid,Vq',cscale)
        hold on;
scatter(new_shape_params(i).old_locations(:,2),new_shape_params(i).old_locations(:,1) ,'MarkerEdgeColor',color_list(i+1,:)) % stimulated locations
  ylabel('x');    xlabel('y');title(['True shape; z = ' num2str(this_z)]) 
colorbar('eastoutside')

%% Posterior variances
figure(2)

tmp=[];
for i=1:length(param_list)
    tmp=[post_var_shape{i};tmp];
end

cscale=[min(tmp) max(tmp)];

 color_list= flip(copper(length(param_list)+1));
 
for i=1:length(param_list)
 subplot(1,length(param_list),i) 
Vq_fits = griddata(new_shape_params(i).locations(:,1),new_shape_params(i).locations(:,2),post_var_shape{i},Xq,Yq);
imagesc(y_grid,x_grid,Vq_fits',cscale)
hold on;
%     caxis(cscale)
% hold on;
if i==1
    ylabel('x');xlabel('y');title(['Prior variance'])  
else
    scatter(new_shape_params(i).old_locations(:,2),new_shape_params(i).old_locations(:,1) ,'MarkerEdgeColor','r') % stimulated locations

 ylabel('x');xlabel('y');title(['Batch ' num2str(i-1)])  
end
colorbar('eastoutside')
end
%% Contrasts
 % Prior to truth:
 figure(3)
 cscale=   [-0.2 0.2];
 color_list= flip(copper(length(param_list)+1));

Vqprime=griddata(new_shape_params(1).locations(:,1),new_shape_params(1).locations(:,2),post_mean_shape{1},Xq,Yq);

Vq = griddata(xy_locations(:,1),xy_locations(:,2),true_shape_values,Xq,Yq);

 for i=1:length(param_list)
 subplot(1,length(param_list),i) 
Vq_fits = griddata(new_shape_params(i).locations(:,1),new_shape_params(i).locations(:,2),post_mean_shape{i},Xq,Yq);
imagesc(y_grid,x_grid,(Vq_fits-Vq)',cscale)
hold on;
%     caxis(cscale)
% hold on;

if i==1
    diff_mat=Vq_fits-Vq;
    diff_vec=[diff_mat(:)];
    diff_vec=diff_vec(~isnan(diff_vec));
    prior_deviation=mean( (diff_vec).^2);
    ylabel('x');xlabel('y');title(['Prior shape - true shape '])  
    
else
    scatter(new_shape_params(i).old_locations(:,2),new_shape_params(i).old_locations(:,1) ,'MarkerEdgeColor','r') % stimulated locations

    diff_mat=Vq_fits-Vq;
    diff_vec=[diff_mat(:)];
    diff_vec=diff_vec(~isnan(diff_vec));
    mean_deviation=mean( (diff_vec).^2);
reduct=100*(1- (mean_deviation/prior_deviation));
 ylabel('x');xlabel('y');title(['Batch ' num2str(i-1) ' Dev. Reduction ' num2str(round(reduct,0)) '%'])  
end

colorbar('eastoutside')
end

%  subplot(1,length(param_list),i) 
% imagesc(y_grid,x_grid,(Vq-Vqprime)',cscale)
%         hold on;
% scatter(new_shape_params(i).old_locations(:,2),new_shape_params(i).old_locations(:,1) ,'MarkerEdgeColor','r') % stimulated locations
%   ylabel('x');    xlabel('y');title(['True shape; z = ' num2str(this_z)]) 
      
% colorbar('eastoutside')
%% Draw the estimated shape values v.s. the true shape values:
% shape_sample_mean
% post_mean=  new_shape_params(i).post_mean(1:size(new_shape_params(i).old_locations,1));

post_mean=cell([length(param_list) 1]);
shape_truth=cell([length(param_list) 1]);
for i = 1:length(param_list)
    post_mean{i}=post_mean_shape{i}(1:size(new_shape_params(i).old_locations,1));
    % post_mean(:,i)=new_shape_params(i).post_mean(1:size(new_shape_params(i).old_locations,1));
    
    shape_truth{i}=zeros(size(new_shape_params(i).old_locations,1),1);
    for i_loc = 1:size(new_shape_params(i).old_locations,1)
        rel_loc=new_shape_params(i).old_locations(i_loc,:);
        this_size = griddata(params.mesh_grid(:,1),params.mesh_grid(:,2),params.mesh_grid(:,3),...
            one_neuron.truth.shape,rel_loc(1),rel_loc(2),rel_loc(3),'linear');
        shape_truth{i}(i_loc)=this_size;
    end
    
end
figure(4)
subplot(1,3,1)
for i=1:length(param_list)
scatter(shape_truth{i},  post_mean{i},'MarkerEdgeColor',color_list(i,:),'MarkerFaceColor', color_list(i,:))
hold on;
line([0 1], [0 1])
xlabel('True shapes')
ylabel('Prior (red) and Posterior (blue) mean)')
xlim([0 1]); ylim([0 1]);
% Vq = griddata(xy_locations(:,1),xy_locations(:,2),true_shape_values,Xq,Yq);
end 
subplot(1,3,2)
for i=1:length(param_list)
scatter(shape_truth{i},abs(post_mean{i}-shape_truth{i}),'MarkerEdgeColor',color_list(i,:),'MarkerFaceColor', color_list(i,:))
hold on;
line([0 1], [0 0])
xlabel('True shapes')
ylabel('Absoluate difference')
% Vq = griddata(xy_locations(:,1),xy_locations(:,2),true_shape_values,Xq,Yq);
end 
subplot(1,3,3)
for i=1:length(param_list)
scatter(post_var_shape{i}(1:size(new_shape_params(i).old_locations,1)), abs(post_mean{i}-shape_truth{i}),...
    'MarkerEdgeColor',color_list(i,:),'MarkerFaceColor', color_list(i,:))
hold on;
% line([0 1], [0 0])
xlabel('Posterior shape variance')
ylabel('Absoluate difference')
% Vq = griddata(xy_locations(:,1),xy_locations(:,2),true_shape_values,Xq,Yq);
end 


end