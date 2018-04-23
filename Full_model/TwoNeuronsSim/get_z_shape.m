function  [z_shape_function, z_dictionary,z_variance]=get_z_shape(path)
load(path);
% Remove the first one because of the singularity 
% Rescale the current:
for i = 2:length(result_z)
    result_z(i).scaled_curr=result_z(i).max_curr/max(result_z(i).max_curr);
end

z_grid = unique([result_z(:).current_z_pos]);
z_shape = z_grid;z_shape(:)=0;
z_shape_list=cell(length(z_grid),1);
z_demeaned_list=cell(length(z_grid),1);

z_grid_list=cell(length(z_grid),1);
z_dictionary=cell(length(result_z)-1,1);
z_count=z_grid;z_count(:)=0;

for i = 2:length(result_z)
    z_dictionary{i-1}=fit(result_z(i).current_z_pos',result_z(i).scaled_curr,'smoothingspline','SmoothingParam',0.07);
    
    true_shape_func =  z_dictionary{i-1};
    true_grid=result_z(i).current_z_pos;
    smoothed_current=true_shape_func(true_grid);
    demeaned_current=result_z(i).scaled_curr-smoothed_current;
    
    for l=1:length(result_z(i).scaled_curr)
        grid_index=find(z_grid==result_z(i).current_z_pos(l));
        z_shape_list{grid_index}=[z_shape_list{grid_index} result_z(i).scaled_curr(l)]; 
        z_demeaned_list{grid_index}=[z_demeaned_list{grid_index} demeaned_current(l)];
        z_grid_list{grid_index}=[z_grid_list{grid_index} z_grid(grid_index)];
        z_shape(grid_index)=z_shape(grid_index)+result_z(i).scaled_curr(l);
        z_count(grid_index)=z_count(grid_index)+1;
    end
end
z_shape_function = fit([z_grid_list{:}]',[z_shape_list{:}]','smoothingspline','SmoothingParam',0.07);
%% Estimate the variance:

sq_dev=([z_shape_list{:}]'-z_shape_function([z_grid_list{:}]')).^2;
z_variance = fit([z_grid_list{:}]',sq_dev,'smoothingspline','SmoothingParam',0.07);


