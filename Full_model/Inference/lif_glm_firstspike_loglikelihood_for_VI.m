function [loss] = lif_glm_firstspike_loglikelihood_for_VI(...
    mpp_this_trial,...
                prob_this_trial,background_rate)
            
%             mpp_this_trial=mpp(i_trial);
             
            n_grid=size(prob_this_trial,2);
        if isempty(mpp_this_trial.times)==true
            temp_prob=...
                bsxfun(@minus,1,sum(prob_this_trial,2));
            singular_index=find(temp_prob< 1e-8); % avoid singularity
            temp_prob(singular_index)=1e-8;
            loss= -n_grid*background_rate+sum(log(temp_prob));
        else
            temp = bsxfun(@plus,background_rate,sum(prob_this_trial(:,round(mpp_this_trial.times)),1));
            
            loss= sum(log(temp));
        end
    end
    
    
    
    