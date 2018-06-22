function [values] = quick_match(new_grid,params)
values=zeros(length(new_grid),1);
for i = 1:length(new_grid)
    [~,im]=min( abs(params.grid -new_grid(i)));
    values(i)=params.values(im);
end

