function plot_max_proj_w_locs(image,locations,experiment_setup)

locs_img = locations(:,1:3);
locs_img(:,1:2) = locs_img(:,1:2)/experiment_setup.image_um_per_px;
locs_img(:,3) = locs_img(:,3)/experiment_setup.stack_um_per_slice;
locs_img = bsxfun(@plus,locs_img,[experiment_setup.image_zero_order_coord' 0]);%
locs_img(:,1:2) = locs_img(:,[2 1]);

% sizes = 5*ones(size(locs_img,1),1);

imagesc(max(image,[],3));
% caxis([0 4000])
M = 12500;
G = linspace(0,1,M)';
myGmap = horzcat(G, zeros(size(G)) , zeros(size(G)));
colormap(myGmap)
axis image;
axis off
hold on;

scatter(locs_img(:,1),locs_img(:,2),5,'w','filled')


    
    