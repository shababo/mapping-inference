function colmat = get_colmat_from_structarray(structarray,fieldname)

colmat = reshape([structarray.(fieldname)],size(structarray(1).(fieldname),1),[])';