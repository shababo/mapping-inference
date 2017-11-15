function rowmat = get_rowmat_from_structarray(structarray,fieldname)

rowmat = reshape([these_trials.cell_IDs],length(structarray.(fieldname)),[])';