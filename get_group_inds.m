function [inds, cell_IDs] = get_group_inds(neighbourhood,group_ID,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    b_id=varargin{1};
else
    b_id=length(neighbourhood.neurons(1).group_ID);
end


inds = arrayfun(@(x) strcmp(x.group_ID{b_id},group_ID), neighbourhood.neurons);
all_cell_IDs = [neighbourhood.neurons.cell_ID];
cell_IDs = all_cell_IDs(inds);