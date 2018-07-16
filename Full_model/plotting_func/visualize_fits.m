function [figure_handle]=visualize_fits(figure_handle,neighbourhood,experiment_setup,varargin)
%%

neurons=neighbourhood.neurons;
n_cell = length(neurons);
flds=fieldnames(neurons(1).posterior_stat(end));
i_batch = length(neurons(1).posterior_stat);
contrast=struct;
simulation_params=experiment_setup.sim;
truth_names = {'optical_gain' 'PR' 'delay_mean' 'delay_var' 'shape'};
for i_field = 1:length(flds)
    this_field = flds{i_field};truth_field = truth_names{i_field};
    if ~strcmp(this_field, 'shapes')
        properties={this_field};summary_stat={'mean'};
        temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);
        contrast.(this_field)=zeros(length(temp_output.(this_field).mean),2);
        contrast.(this_field)(:,1)=temp_output.(this_field).mean;
        for i_cell = 1:length(neurons)
            contrast.(this_field)(i_cell,2)= neurons(i_cell).truth.(truth_field);
        end
        if strcmp(truth_field, 'delay_var')
            contrast.(this_field)(i_cell,2)=sqrt(contrast.(this_field)(i_cell,2));
        end
    else % if this is the shape:
        contrast.(this_field).values=cell([n_cell 1]);
        for i_cell = 1:n_cell
            if neurons(i_cell).truth.PR>0
                est_shapes = neurons(i_cell).posterior_stat(i_batch).shapes;
                true_shape=neurons(i_cell).truth.shape;
                this_locs=est_shapes.locations;
                contrast.(this_field).values{i_cell}=zeros(length(est_shapes.mean),2);
                contrast.(this_field).values{i_cell}(:,1)=est_shapes.mean;
                this_size = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
                    true_shape,this_locs(:,1),this_locs(:,2),this_locs(:,3));
                this_size(isnan(this_size))=0;
                contrast.(this_field).values{i_cell}(:,2)=this_size;
            end
        end
    end
end

%% plotting:
figure(figure_handle);
% figure(1)
connected_cells = contrast.PR(:,2)>0;
color_set={'r' 'b'};
color_for_shapes= lines(sum(connected_cells));
i_count = 1;
fig_index = [1 3 4 2];
for i_field = 2:length(flds)
    this_field = flds{i_field};
    h(i_field - 1)=subplot(2,2,fig_index(i_field - 1));
 
%     figure(i_field)
       
    if ~strcmp(this_field, 'shapes')
        
        scatter(contrast.(this_field)(connected_cells,1),contrast.(this_field)(connected_cells,2),...
            'MarkerFaceColor',color_set{2},'MarkerEdgeColor',color_set{2});
        hold on;
        scatter(contrast.(this_field)(~connected_cells,1),contrast.(this_field)(~connected_cells,2),...
              'MarkerFaceColor',color_set{1},'MarkerEdgeColor',color_set{1});
        hold on;
        
        for i_cell = 1:n_cell
            text(contrast.(this_field)(i_cell,1),contrast.(this_field)(i_cell,2),num2str(i_cell))
        end
       xmax=max([contrast.(this_field)(:,1)]);
       ymax=max([contrast.(this_field)(:,2)]);
       switch this_field 
           case 'delay_mu'
               title_txt = 'delay mean';
           case 'delay_sigma'
               title_txt = 'delay sigma';
           case 'PR'
               title_txt = [this_field ', Batch ' num2str(i_batch-1)]; 
           otherwise
               title_txt=this_field;
       end
       title(title_txt);
    else
        xmax=0;ymax=0;
        for i_cell = 1:n_cell
            if neurons(i_cell).truth.PR>0
                temp_shape=contrast.(this_field).values{i_cell};
                gains= contrast.gain(i_cell,:);
                scatter(temp_shape(:,1)*gains(1),temp_shape(:,2)*gains(2),...
                    'MarkerFaceColor',color_for_shapes(i_count,:),'MarkerEdgeColor',color_for_shapes(i_count,:));
                hold on;
                i_count = i_count+1;
                xmax=max([xmax; temp_shape(:,1)*gains(1)]);
                ymax=max([ymax; temp_shape(:,2)*gains(2)]);
            text(temp_shape(:,1)*gains(1),temp_shape(:,2)*gains(2),num2str(i_cell))
        
            end
            
        end
        title(['Gain * shape'])
    end
    
    plot([0 100], [0 100],'Color','k')
    hold off;
    xlim([0 xmax])
    ylim([0 ymax])
    xlabel('Estimates')
    xlabel('True values')
end



