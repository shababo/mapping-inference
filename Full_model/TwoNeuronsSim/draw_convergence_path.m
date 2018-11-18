function [parameter_path]=draw_convergence_path(parameter_history,elbo_rec, params)

% params.truth store the true values for neurons
% params.plot.do a binary for whether drawing the plots
% params.plot.handles handles for figures.


% params =struct;
% params.truth=struct;
% params.truth.neurons=neurons;
% params.plot.do=true;
% params.plot.original_scale=true; % whether to plot values in the original scales 
%% Record the path of fits
% This is basically a reformat of the parameter_history
parameter_path =struct;
n_iterations=size(parameter_history,1);
n_neurons=size(parameter_history,2);
fldnames = fieldnames(parameter_history(1,1));
for i=1:length(fldnames)
    i_field = fldnames{i};
    if ~strcmp(i_field,'shapes')
        tmp_mat = zeros(n_iterations,n_neurons); tmp_mat_sig= tmp_mat;
        for i_neuron = 1:n_neurons
            for i_iteration=1:n_iterations
                tmp_mat(i_iteration,i_neuron) = parameter_history(i_iteration,i_neuron).(i_field).mean;
                tmp_mat_sig(i_iteration,i_neuron) = exp(parameter_history(i_iteration,i_neuron).(i_field).log_sigma);
            end
        end
        parameter_path.(i_field).mean=tmp_mat; parameter_path.(i_field).sigma=tmp_mat_sig;
        parameter_path.(i_field).bounds=parameter_history(i_iteration,i_neuron).(i_field).bounds;
    else
        parameter_path.(i_field)=struct;
        parameter_path.(i_field).neurons(n_neurons)=struct;
        for i_neuron = 1:n_neurons
            tmp_mat = zeros(n_iterations,length(parameter_history(1,i_neuron).(i_field).mean));
            tmp_mat_sig= tmp_mat;
            for i_iteration=1:n_iterations
                tmp_mat(i_iteration,:) = parameter_history(i_iteration,i_neuron).(i_field).mean;
                tmp_mat_sig(i_iteration,:) = exp(parameter_history(i_iteration,i_neuron).(i_field).log_sigma);
            end
            parameter_path.(i_field).neurons(i_neuron).mean=tmp_mat;
            parameter_path.(i_field).neurons(i_neuron).sigma=tmp_mat_sig;
            parameter_path.(i_field).neurons(i_neuron).location=parameter_history(1,i_neuron).shapes.locations;
            parameter_path.(i_field).neurons(i_neuron).bounds=parameter_history(1,i_neuron).shapes.bounds;
        end
    end
end
parameter_path.elbo_path = elbo_rec;
%% Record the true values
draw_truth=false;
if isfield(params,'truth')
    draw_truth=true;
end

%% Draw plots
if params.plot.do
    % Figure 1: Convergence path for all neurons for all pathes except for the shapes
    figure(3)
    for i=1:length(fldnames)
        i_field = fldnames{i};
        if ~strcmp(i_field,'shapes')
            for i_neuron =1 :n_neurons
                i_plot = i_neuron+(i-1)*n_neurons;
                subplot(length(fldnames)-1,n_neurons,i_plot);
                
                % Should wrap this up as a function:
                if params.plot.original_scale
                    means_tmp=parameter_path.(i_field).mean(:,i_neuron);
                    sigmas_tmp=parameter_path.(i_field).sigma(:,i_neuron);
                    means = exp(means_tmp)./(1+exp(means_tmp))*...
                        (parameter_path.(i_field).bounds.up-parameter_path.(i_field).bounds.low)+parameter_path.(i_field).bounds.low;
                    mps = exp(means_tmp+sigmas_tmp)./(1+exp(means_tmp+sigmas_tmp))*...
                        (parameter_path.(i_field).bounds.up-parameter_path.(i_field).bounds.low)+parameter_path.(i_field).bounds.low;
                    mms = exp(means_tmp-sigmas_tmp)./(1+exp(means_tmp-sigmas_tmp))*...
                        (parameter_path.(i_field).bounds.up-parameter_path.(i_field).bounds.low)+parameter_path.(i_field).bounds.low;
                       if draw_truth
                    true_value = params.truth.neurons(i_neuron).truth.(i_field)*ones(n_iterations,1);
                       end
                else
                    means_tmp=parameter_path.(i_field).mean(:,i_neuron);
                    sigmas_tmp=parameter_path.(i_field).sigma(:,i_neuron);
                    means =mean_tmp;
                    mps = means_tmp+sigmas_tmp;
                    mms = means_tmp-sigmas_tmp;
                     if draw_truth
                    
                    true_tmp=params.truth.neurons(i_neuron).truth.(i_field);
                    true_value = log(true_tmp./(1-true_tmp))*ones(n_iterations,1);
                     end
                end
                %------------------------------%
                
                
                plot(1:n_iterations, means,'color','k')
                hold on;
                line(1:n_iterations, mps,'LineStyle',':','color','k')
                hold on;
                line(1:n_iterations, mms,'LineStyle',':','color','k')
                hold on;
                title_string= ['Neuron ' num2str(i_neuron)];
                if draw_truth
                    line(1:n_iterations, true_value,'color','r','LineStyle','--')
                   title_string=['Neuron ' num2str(i_neuron) '; PR = ' num2str(params.truth.neurons(i_neuron).truth.PR ) ];
                end
                xlabel('Iteration')
                ylabel(i_field)
                
                title(title_string)
            end
        else
            % Draw shapes:
% %             Use dotted curves to represent paths of spikes
%             for i_neuron =1 :n_neurons
%                 i_plot = i_neuron+(i-1)*n_neurons;
%                 subplot(length(fldnames),n_neurons,i_plot);
%                 n_shapes = size(parameter_path.(i_field).neurons(i_neuron).mean,2);
%                 color_list= lines(n_shapes);
%                 for i_shape = 1:n_shapes
%                     
%                      % Should wrap this up as a function:
%                 if params.plot.original_scale
%                     means_tmp=parameter_path.(i_field).neurons(i_neuron).mean(:,i_shape);
%                     sigmas_tmp=parameter_path.(i_field).neurons(i_neuron).sigma(:,i_shape);
%                     ub =parameter_path.(i_field).neurons(i_neuron).bounds.up(i_shape);
%                     lb =parameter_path.(i_field).neurons(i_neuron).bounds.low(i_shape);
%                     
%                     means = exp(means_tmp)./(1+exp(means_tmp))*(ub-lb)+lb;
%                     mps = exp(means_tmp+sigmas_tmp)./(1+exp(means_tmp+sigmas_tmp))*(ub-lb)+lb;
%                     mms = exp(means_tmp-sigmas_tmp)./(1+exp(means_tmp-sigmas_tmp))*(ub-lb)+lb;
%  
%                     rel_loc =  parameter_path.(i_field).neurons(i_neuron).location(i_shape,:);
%                     this_size = griddata(params.mesh_grid(:,1),params.mesh_grid(:,2),params.mesh_grid(:,3),...
%                         params.truth.neurons(i_neuron).truth.shape,rel_loc(1),rel_loc(2),rel_loc(3),'linear');
%                     true_value = this_size*ones(n_iterations,1);
%                 else
%                     means_tmp=parameter_path.(i_field).mean(:,i_neuron);
%                     sigmas_tmp=parameter_path.(i_field).sigma(:,i_neuron);
%                     means =mean_tmp;
%                     mps = means_tmp+sigmas_tmp;
%                     mms = means_tmp-sigmas_tmp;
%                     
%                     ub =parameter_path.(i_field).neurons(i_neuron).bounds.up(i_shape);
%                     lb =parameter_path.(i_field).neurons(i_neuron).bounds.low(i_shape);
%                     rel_loc =  parameter_path.(i_field).neurons(i_neuron).location(i_shape,:);
%                     this_size = griddata(params.mesh_grid(:,1),params.mesh_grid(:,2),params.mesh_grid(:,3),...
%                         params.truth.neurons(i_neuron).truth.shape,rel_loc(1),rel_loc(2),rel_loc(3),'linear');
%                     
%                     
%                     true_value = log(this_size./(1-this_size))*ones(n_iterations,1);
%                 end
%                 %------------------------------%
%                 
%                     if draw_truth
%                         
%                     plot([true_value(1) true_value(1)], means([1 end]),'color',color_list(i_shape,:),...
%                         'LineStyle','--')
%                     hold on;
%                     scatter(true_value(1), means(end),'MarkerFaceColor',color_list(i_shape,:))
%                     hold on;
%                   
% %                     line(1:n_iterations, mps,'color',color_list(i_shape,:),'LineStyle',':')
% %                     hold on;
% %                     line(1:n_iterations, mms,'color',color_list(i_shape,:),'LineStyle',':')
% %                     hold on;
%                     title_string = ['Neuron ' num2str(i_neuron)];
%                     
% %                     line(1:n_iterations, true_value,'LineStyle','--','color',color_list(i_shape,:))
% %                         hold on;
% %                         title_string=['Neuron ' num2str(i_neuron) '; PR = ' num2str(params.truth.neurons(i_neuron).truth.PR ) ];
%                     end
%                 end
%                     xlabel('Iteration')
%                 ylabel(i_field)
%            
%              title(title_string)
%                line([0 1],[0, 1])
%             end
        end
    end
    
    % Figure 2: a scatter plot between ELBO and some parameters
    figure(8)
    for i=1:length(fldnames)
        i_field = fldnames{i};
        if ~strcmp(i_field,'shapes')
            for i_neuron =1 :n_neurons
                i_plot = 2*i_neuron-1+(i-1)*2*n_neurons;
                subplot(length(fldnames)-1,2*n_neurons,i_plot);
                sigmas=parameter_path.(i_field).sigma(:,i_neuron);
                
                     % Should wrap this up as a function:
                if params.plot.original_scale
                      means_tmp=parameter_path.(i_field).mean(:,i_neuron);
                    ub =parameter_path.(i_field).bounds.up;
                    lb =parameter_path.(i_field).bounds.low;
                    means = exp(means_tmp)./(1+exp(means_tmp))*(ub-lb)+lb;
                    if draw_truth
                    true_value = params.truth.neurons(i_neuron).truth.(i_field)*ones(2,1);
                    end
                    else
                           means_tmp=parameter_path.(i_field).mean(:,i_neuron);
                    ub =parameter_path.(i_field).bounds.up;
                    lb =parameter_path.(i_field).bounds.low;
                  means =mean_tmp;
                  if draw_truth
                        true_value = log(true_tmp./(1-true_tmp))*ones(2,1);
                  end  
              end
                %------------------------------%
            
                
                plot(means,parameter_path.elbo_path(2:end))
                hold on;
                
                  title_string = ['Neuron ' num2str(i_neuron)];
                  
                if draw_truth
                    line(true_value,[min(parameter_path.elbo_path), max(parameter_path.elbo_path)],'color','r');
                              title_string=['Neuron ' num2str(i_neuron) '; PR = ' num2str(params.truth.neurons(i_neuron).truth.PR ) ];
                hold on;
                end
                xlabel([i_field ' (mean)'])
                ylabel('ELBO')
                title(title_string)
                i_plot = 2*i_neuron+(i-1)*2*n_neurons;
                subplot(length(fldnames)-1,2*n_neurons,i_plot);
                plot(sigmas,parameter_path.elbo_path(2:end))
                title(title_string)
                
                xlabel([i_field ' (sigma)'])
                ylabel('ELBO')
                
            end
        end
        
    end
end
