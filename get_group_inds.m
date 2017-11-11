function inds = get_group_inds(neighbourhood,group_ID)

inds = arrayfun(@(x) strcmp(x.group_ID,group_ID), neighbourhood.neurons);