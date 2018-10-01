
function [figure_handle]=visualize_fits_batch(figure_handle,neurons,simulation_params)
%%

% neurons=neighbourhood.neurons;
n_cell = length(neurons);
% flds=fieldnames(neurons(1).posterior_stat(end));
flds={'gain' 'shapes'};
i_batch = length(neurons(1).posterior_stat);
contrast=struct;
truth_names = {'optical_gain' 'shape'};
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
                
                if ~isempty(est_shapes.mean)
                    
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
end

%% plotting:
figure(figure_handle);
% figure(1)
color_set={'r'};
color_for_shapes= copper(1);
i_count = 1;

 this_field = 'shapes';
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
    
    plot([0 100], [0 100],'Color','k')
    hold off;
    xlim([0 xmax])
    ylim([0 ymax])
    xlabel('Estimates')
    ylabel('True values')