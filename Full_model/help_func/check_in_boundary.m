function [relevant_flag]= check_in_boundary(rel_position,boundary_params)
%%
% boundary_params = [30 30 70];
relevant_flag=true;

scaled_distance=rel_position./boundary_params;
relevant_flag = (sum(scaled_distance.^2 )< 1);
