function cell_features = filter_nuclear_detection(cell_features, size_lims)

bad_ids = zeros(size(cell_features,2),1);
for i = 1:size(cell_features,2)
    if cell_features(1,i) < 35%cell_features(5,i) > size_lims || cell_features(8,i) > size_lims || ...
        %cell_features(5,i) < 2.5 || cell_features(8,i) < 2.5  || ...
        %any(cell_features([2 3 4],i) < 0) || ...
        
        
        bad_ids(i) = 1;
    end
end

cell_features(:,logical(bad_ids)) = [];