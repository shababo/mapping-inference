function [ cluster_of_cells] = find_clusters(...
    stim_size, cell_list,stim_threshold)
        % Break the trials into unrelated clusters
        stim_size=stim_size(:,cell_list);
        if length(cell_list)>1
            % Use inner product:
            stim_size=(stim_size>stim_threshold).*stim_size;
            adj_corr= abs(stim_size'*stim_size)./size(stim_size,1);
%             adj_corr=1*(adj_corr> stim_threshold^2);
            adj_corr=adj_corr>0;
            adj_corr=adj_corr+diag(ones(size(adj_corr,1),1));
            
            cc_corr=expm(adj_corr);
            cell_cluster_ind=zeros(length(cell_list),1);
            cluster_id=1;
            for i_cell_idx = 1:length(cell_list)
                if cell_cluster_ind(i_cell_idx)==0
                    this_id=cluster_id;
                    cluster_id=cluster_id+1;
                else
                    this_id = cell_cluster_ind(i_cell_idx);
                end
                cell_cluster_ind( find(cc_corr(:,i_cell_idx)))=this_id;
            end
            % Now turn the cell_cluster_ind into list of cells
            n_cluster=max(cell_cluster_ind);
            cluster_of_cells= cell([n_cluster 1]);
            for i_cluster = 1:n_cluster
                cluster_of_cells{i_cluster}=find(cell_cluster_ind==i_cluster);
            end
            
        else
            cluster_of_cells= cell([1 1]);
            cluster_of_cells{1}=1;
        end
        
      
        
end




