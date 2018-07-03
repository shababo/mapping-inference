function [output]=grab_values_from_neurons(i_batch,neurons,properties,summary_stat)

number_of_cells=length(neurons);
output=struct;
for i_properties = 1:length(properties)
    this_property=properties{i_properties};
    output.(this_property)=struct;
    for i_stat = 1:length(summary_stat)
        this_stat=summary_stat{i_stat};
        temp_values=zeros(number_of_cells,1);
        for i_cell = 1:number_of_cells
            
            % add a statement to check if we have stop updating this cell 
            if i_batch > length(neurons(i_cell).posterior_stat)
                temp_values(i_cell)= neurons(i_cell).posterior_stat(end).(this_property).(this_stat);    
            else
                temp_values(i_cell)= neurons(i_cell).posterior_stat(i_batch).(this_property).(this_stat);    
            end
        end
        output.(this_property).(this_stat)=temp_values;
    end
end
end
