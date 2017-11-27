function rowmat = get_rowmat_from_structarray(structarray,fieldname)

rowmat = reshape([structarray.(fieldname)],size(structarray(1).(fieldname),2),[])';