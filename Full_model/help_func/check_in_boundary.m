function [relevant_flag]= check_in_boundary(rel_position,boundary_params)
%%
% boundary_params = [30 30 70];
relevant_flag=true;
for i_coord = 1:length(rel_position)
    relevant_flag = relevant_flag & (abs(rel_position(i_coord))<boundary_params(i_coord));
end
