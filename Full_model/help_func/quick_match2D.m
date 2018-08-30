function [values] = quick_match2D(new_grid,params)
values=zeros(size(new_grid,1),1);
for i = 1:size(new_grid,1)
    [~,im]=min( sum(   [(params.grid -new_grid(i,:)).^2]'  ));
    values(i)=params.values(im);
end

