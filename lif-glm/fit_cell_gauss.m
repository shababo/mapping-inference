function [sigma] = fit_cell_gauss(cell_shape,locations,sig0)
sigma = zeros(3,1);
for i = 1:length(sigma)
    f = fit(x.',y.','gauss2')
    plot(f,x,y)

% cell_shape = gauss3D_fit(sig0,locations);
% assignin('base','cell_shape_test',cell_shape)
% % return
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare';
% sigma = nlinfit(locations,cell_shape(:),@gauss3D_fit,sig0(:),opts);
% sigma = reshape(sigma,3,3);


end

function cell_shape = gauss3D_fit(sigma,locations)
z_depths = sort(unique(locations(:,3)));
cell_shape = zeros(9*9*7,1);
for i = 1:size(locations,1)
%     locations(i,:)
    inds = locations(i,[1 2])/10+5; %HARD CODED FOR NOW :(
    inds(3) = find(locations(i,3) == z_depths);
    cell_shape(sub2ind([9 9 7],inds(1),inds(2),inds(3))) = ...
        exp(-0.5*(locations(i,:)*pinv(reshape(sigma,3,3))*locations(i,:)'));
%         mvnpdf(locations(i,:),[0 0 0],reshape(sigma,3,3))/mvnpdf([0 0 0],[0 0 0],reshape(sigma,3,3));
end
end



lininds = sub2ind([9 9 7],sorted_shape_locs(:,1)/10+5,sorted_shape_locs(:,2)/10+5,arrayfun(@(x) find(x == z_depths),sorted_shape_locs(:,3)));