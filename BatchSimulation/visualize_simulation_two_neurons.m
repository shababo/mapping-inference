function [ ]=visualize_simulation_two_neurons(fig_handle,result_mat,list_one,list_two,i_conn,axis_choice,native_indicator)
%%
n_grid = 400;
figure(fig_handle)
mat_one = zeros(length(list_one),length(list_two));
for i_gain = 1:length(list_one)
    for i_dist = 1:length(list_two)
        mat_one(i_gain,i_dist)=mean(result_mat(i_gain,i_conn,axis_choice,i_dist,:));
    end
end

grid_one=min(list_one):range(list_one)/n_grid:max(list_one);
grid_two=min(list_two):range(list_two)/n_grid:max(list_two);
[mesh_one,mesh_two]=meshgrid(list_one,list_two);
[Oneq,Twoq] = meshgrid(grid_one,grid_two);
vec_one=reshape(mesh_one,[],1);vec_two=reshape(mesh_two,[],1);
mat_vec=reshape(mat_one,[],1);
for i_vec = 1:length(mat_vec)
    this_one=vec_one(i_vec);this_two=vec_two(i_vec);
    mat_vec(i_vec)=mat_one(find(list_one==this_one), find(list_two==this_two));
end

Vq = griddata(vec_one,vec_two,flip(mat_vec,2),Oneq,Twoq);

% Vq = griddata(list_one,list_two',mat_one,Oneq,Twoq);
%
cscale=[0 1];
if native_indicator
imagesc(grid_one,grid_two,flip(mat_one,2)',cscale)
else
imagesc(grid_one,grid_two,Vq,cscale)
end
colorbar()

