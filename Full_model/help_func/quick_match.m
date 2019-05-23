function [values] = quick_match(new_grid,params)
values=zeros(size(new_grid,1),1);
for i = 1:length(new_grid)
if size(new_grid,2)>1
[~,im]=min( sum( abs(params.grid -new_grid(i,:)).^2,2) );
else
    [~,im]=min(abs(params.grid -new_grid(i)) );
end
    values(i)=params.values(im);
end

