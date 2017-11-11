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
            temp_struct=neurons(i_cell).(this_property);
            temp_values(i_cell)= temp_struct(i_batch).(this_stat);
        end
        output.(this_property).(this_stat)=temp_values;
    end
end
end

end