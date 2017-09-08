function cell_features = filter_nuclear_detection(cell_features, size_lims)

% bad_ids = zeros(size(cell_features,2),1);
% for i = 1:size(cell_features,2)
%     if cell_features(1,i) < 20 ...
%         ...%|| cell_features(size_min,i) > size_lims || cell_features(8,i) > size_lims ...
%         ...%|| cell_features(size_min,i) < size_min || cell_features(8,i) < size_min ...
%         || any(cell_features([2 3 4],i) < 0) ...
%         || cell_features(4,i) < 2
%         
%         
%         bad_ids(i) = 1;
%     end
% end
% 
% cell_features(:,logical(bad_ids)) = [];

size_max = 40;
size_min = 0;
% cell_features(:,cell_features(5,:)>size_max) = [];
% cell_features(:,cell_features(5,:)<size_min) = [];
% cell_features(:,cell_features(8,:)<size_min) = [];
% cell_features(:,cell_features(10,:)>size_max) = [];
% cell_features(:,cell_features(10,:)<size_min) = [];
cell_features(:,cell_features(1,:)<140) = []; %140
