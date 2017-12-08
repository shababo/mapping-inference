function [inds, cell_IDs] = get_group_inds(neighbourhood,group_ID)
b_id=neighbourhood.batch_ID;
inds = arrayfun(@(x) strcmp(x.group_ID{b_id},group_ID), neighbourhood.neurons);
all_cell_IDs = [neighbourhood.neurons.cell_ID];
cell_IDs = all_cell_IDs(inds);