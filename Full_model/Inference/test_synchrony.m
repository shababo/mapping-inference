function [disconnected_indicators] = test_synchrony(disconnected_list, connected_list,...
    cell_list_map,mpp_connected,mpp_undefined,loc_to_cell,loc_to_cell_nuclei,background_rate,gamma_bound)
%  alpha=parameter_history.alpha(:,end);
%  beta=parameter_history.beta(:,end);


%     disconnected_list=disconnected_cells{iter};
%     connected_list=alive_cells{iter}+connected_cells{iter};
%             
            cell_list= find(disconnected_list);
            connected_cells_list=find(connected_list);
        
            disconnected_indicators=ones(length(disconnected_list),1);
        
            for i_cell_index = 1:length(cell_list)
                % For each cell, pick out trials that this cell was
                % targeted
                % then exclude the trials where there are known connected
                % cells involved
                i_cell=cell_list_map(cell_list(i_cell_index));
                mpp_related=struct([]);
                related_index =[];
                for i_trial = 1:length(mpp_connected)
                    if loc_to_cell_nuclei(mpp_connected(i_trial).locations)==i_cell
                        related_index=[related_index i_trial];
                    end
                end
                mpp_related=mpp_connected(related_index);
                related_index =[];
                for i_trial = 1:length(mpp_undefined)
                    this_trial_locations=mpp_undefined(i_trial).locations(~isnan(mpp_undefined(i_trial).locations));
                    temp_indicator=sum( loc_to_cell(this_trial_locations)==i_cell );
                    known_connection=intersect(loc_to_cell(this_trial_locations),...
                        cell_list_map(connected_cells_list));
                    if isempty(known_connection) && temp_indicator
                        related_index=[related_index i_trial];
                    end
                end
                if isempty(mpp_related)
                    mpp_related=mpp_undefined(related_index); 
                else
                    mpp_related(end+(1:length(related_index)))=mpp_undefined(related_index);
                end 
                
                event_times=sort([mpp_related.times]);
                synchrony_threshold = 1/(background_rate*length(mpp_related)); 
                % mean time gap if these are background events 
                if length(event_times)> (length(mpp_related)*gamma_bound.low)
                    time_gaps = event_times(2:end)-event_times(1:(end-1));
                   if mean(time_gaps) < synchrony_threshold/2 % a simple threshold.. need to replace with testing
                       disconnected_indicators(cell_list(i_cell_index))=0; % not disconnected 
                   end
                end
            end
end