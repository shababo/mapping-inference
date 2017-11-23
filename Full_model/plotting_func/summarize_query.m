function [query_summary] = summarize_query(experiment_query)

group_names=setdiff(fieldnames(experiment_query),'batch_ID');

query_summary=struct;
for i_group = 1:length(group_names)
   this_group = group_names{i_group};
   
   
   cell_list=unique([experiment_query.(this_group).trials(:).cell_IDs]);
   mfactor = max(cell_list)+1;
   unique_loc_list=unique([experiment_query.(this_group).trials(:).location_IDs]*mfactor+...
       [experiment_query.(this_group).trials(:).cell_IDs]);
   
   %query_summary.(this_group);
   if isfield(experiment_query.(this_group).trials(1),'event_times')
       count_event=true;
       query_summary(1).(this_group)=zeros(length(unique_loc_list),4);
       % cell ID, loc ID, trial counts, event counts 
   else
       count_event=false;
       query_summary(1).(this_group)=zeros(length(unique_loc_list),3);
       % cell ID, loc ID, trial counts 
   end
   for i_loc = 1:length(unique_loc_list)
      query_summary.(this_group)(i_loc, 1)=mod(unique_loc_list(i_loc),mfactor);
      query_summary.(this_group)(i_loc, 2)=floor(unique_loc_list(i_loc)/mfactor);
   end
   
   for i_trial = 1:length(experiment_query.(this_group).trials)
       loc_temp=experiment_query.(this_group).trials(i_trial).location_IDs*mfactor+...
       experiment_query.(this_group).trials(i_trial).cell_IDs;
       for i_loc= 1:length(loc_temp)
          i_unique_loc = find(unique_loc_list==loc_temp(i_loc) );
          query_summary.(this_group)(i_unique_loc,3)=query_summary.(this_group)(i_unique_loc,3)+1;
          if count_event
              event_ind=length(experiment_query.(this_group).trials(i_trial).event_times)>0;
          query_summary.(this_group)(i_unique_loc,4)=query_summary.(this_group)(i_unique_loc,4)+...
              event_ind;
          end
       end
   end
    
end



end