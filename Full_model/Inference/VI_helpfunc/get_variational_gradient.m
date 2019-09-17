function [this_gradient]=get_variational_gradient(variational_samples,raw_samples,params,varargin)
%%

% if ~isempty(varargin)
%     only_PR = varargin{1};
% else
%     only_PR = false;
% end
% calculate the gradients
fldnames = fieldnames(params(1));
n_cell = length(params);
% clear('this_gradient')
this_gradient(n_cell)= struct;
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
%         if (only_PR & strcmp(fldnames{i_field},'PR')) | (~only_PR)
            
            this_params=params(i_cell).(fldnames{i_field});
            this_sample=variational_samples(i_cell).(fldnames{i_field});
            this_raw_sample = raw_samples(i_cell).(fldnames{i_field});
            if ~strcmp(fldnames{i_field},'shapes') & ~strcmp(fldnames{i_field},'xy') & ~strcmp(fldnames{i_field},'z') 
                switch this_params.dist
                    case {'normal','log-normal', 'logit-normal'}
                        dmean =   (this_raw_sample-this_params.mean)./exp(this_params.log_sigma);
                        dsigma = -1/2+((this_params.mean-this_raw_sample).^2)./(2*exp(this_params.log_sigma));
                    case {'spiked-logit-normal'}
                        dlogit= (this_sample==0)/(1+exp(this_params.prob_logit))-...
                            (this_sample==0)*exp(this_params.prob_logit)/(1+exp(this_params.prob_logit));
                        dmean =   (this_params.mean-this_raw_sample)./exp(2*this_params.log_sigma);
                        dsigma = -1/2+ (this_params.mean-this_raw_sample).^2./exp(2*this_params.log_sigma);
                        this_gradient(i_cell).(fldnames{i_field}).prob_logit=dlogit;
                        
                end
            else % shapes
                switch this_params.dist
                    case {'normal','log-normal', 'logit-normal'}
                        for i_loc = 1:length(this_raw_sample)
                            dmean(i_loc) =   (this_params.mean(i_loc)-this_raw_sample(i_loc))./exp(this_params.log_sigma(i_loc));
                            dsigma(i_loc) = -1/2+ (this_params.mean(i_loc)-this_raw_sample(i_loc)).^2./(2*exp(this_params.log_sigma(i_loc)));
                        end
                    case 'mvn'
                         mean_prod=(this_params.bounds.up-this_params.bounds.low).*exp(this_params.mean)./((1+exp(this_params.mean)).^2);
                         this_mean=(this_params.bounds.up-this_params.bounds.low).*exp(this_params.mean)./(1+exp(this_params.mean))+this_params.bounds.low;
                         dmean=  -mean_prod.*(this_params.Sigma_tilde_inv*(this_mean-this_raw_sample));
                        % dsigma= diag(this_params.Sigma_inv) .*( ((this_mean-this_raw_sample).^2).*exp(-this_params.log_sigma)/2 - diag(this_params.Sigma_tilde).*exp(-this_params.log_sigma)/2);
                        dsigma= ((this_raw_sample-this_mean).^2).*exp(-this_params.log_sigma)/2 - diag(this_params.Sigma_tilde).*exp(-this_params.log_sigma)/2;
                end
            end
            
            this_gradient(i_cell).(fldnames{i_field}).type =  this_params.type;
            this_gradient(i_cell).(fldnames{i_field}).mean = dmean;
            this_gradient(i_cell).(fldnames{i_field}).sigma = dsigma;
            
%         end
    end
end


