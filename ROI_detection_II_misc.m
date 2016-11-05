%% Visualize the map 
figure(10)
R=2500;
x_range = ceil(sqrt(R)/(x_dense(2)-x_dense(1)));
y_range = ceil(sqrt(R)/(y_dense(2)-y_dense(1)));

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*10);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,1);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

    xlim([20,460]);
    ylim([-900,-400]);

%  thresh =0.3;

%     xlim([20,460]);
%     ylim([-900,-400]);
%      potential_neuron_regional = scatter(new_trials(induced_prob>thresh,1),...
%      -new_trials(induced_prob>thresh,2),amplitude_threshold(j+1)*10,...
%     'filled','o');
%     set(potential_neuron_regional,'MarkerFaceColor','r');
%     alpha(potential_neuron_regional,0.5);
%     
% j=2
% threshold_region_vec= zeros(num_dense^2,1);
% for i = 1:num_dense
%     for l = 1:num_dense
%         threshold_region_vec((i-1)*num_dense + l) =likelihood(i,l,j)+1;
%     end
% end
% 
%  thresh =quantile(threshold_region_vec,0.9);
%     xlim([20,460]);
%     ylim([-900,-400]);
%      potential_neuron_regional = scatter(Z_dense(threshold_region_vec>thresh,1),...
%      -Z_dense(threshold_region_vec>thresh,2),amplitude_threshold(j+1)*10,...
%     'filled','o');
%     set(potential_neuron_regional,'MarkerFaceColor','b');
%     alpha(potential_neuron_regional,0.5);
    
for j = 1:(num_threshold-1)       
     potential_neuron_regional = scatter(neurons_regional{j}(:,1),...
     -neurons_regional{j}(:,2),amplitude_threshold(j+1)*10,...
    'filled','o');
    set(potential_neuron_regional,'MarkerFaceColor','r');
    alpha(potential_neuron_regional,0.5);
     
     potential_neuron_local = scatter(neurons_local{j}(:,1),...
         -neurons_local{j}(:,2),amplitude_threshold(j+1)*10,...
         'filled','o');
     set(potential_neuron_local,'MarkerFaceColor','b');
     alpha(potential_neuron_local,0.5);
     hold on
end

hold off
%saveas(10,'../Data/Figure3A.jpg')
view(2)


%% Infer the locations of neurons on the same dense grid
num_dense = 2e2+1;

likelihood = zeros(num_dense,num_dense,num_threshold-1);
x_dense = zeros(num_dense,1);
y_dense = zeros(num_dense,1);
x_dense = (0:(num_dense-1))*(max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1));
y_dense = (0:(num_dense-1))*(max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2));
Z_dense = zeros(num_dense^2,2);
% Take j=5 for example
for j = 1:(num_threshold-1)
    estimates = induced_prob_categories(:,j);
    threshold = 2*induced_prob_categories_sd(:,j);
    % we would need a better standard deviation estimates
    probability = estimates(estimates>threshold);
    prob_sd = induced_prob_categories_sd(estimates>threshold,j);
    locations = new_trials(estimates > threshold,:);
    reference_prob = min(max(probability),1);
    for i = 1:num_dense
        for l = 1:num_dense
            Z_dense((i-1)*num_dense + l,:) = [x_dense(i) y_dense(l)];
            for k = 1:size(probability,1)
                dist_scaled = (locations(k,1)-x_dense(i))^2/A(1,1)+(locations(k,2)- y_dense(l))^2/A(2,2);
                p_ijk = reference_prob*exp(-0.5*dist_scaled);
                % then check the Gaussian density
                likelihood(i,l,j) = likelihood(i,l,j)+normpdf(p_ijk,probability(k),prob_sd(k));
            end
        end
    end
end

%% Find the regional maximum among the dense locations 
% There is a costumized threshold to be chosen...
threshold_level = 3;
quantile_level=0.95;
R=2500;
x_range = ceil(sqrt(R)/(x_dense(2)-x_dense(1)));
y_range = ceil(sqrt(R)/(y_dense(2)-y_dense(1)));

neurons_regional = cell(num_threshold-1,1);
neurons_local = cell(num_threshold-1,1);

for j = 1:num_threshold-1
    likelihood_thresholds =quantile( reshape(likelihood(:,:,j),[num_dense^2 1]), [0.90  0.95 0.99]);
    [cent{j}, varargout(:,:,j)]=FastPeakFind(likelihood(:,:,j),likelihood_thresholds(threshold_level));
    
    neuron_centers = zeros(size(cent{j},1)/2, 2);
    for i = 1:(size(cent{j},1)/2)
        neuron_centers(i,:) = [x_dense(cent{j}( 2*(i-1)+2)) y_dense(cent{j}( 2*(i-1)+1))]; 
    end
    
    threshold_region_vec= zeros(num_dense^2,1);
    for i = 1:num_dense
        for l = 1:num_dense
            threshold_region_vec((i-1)*num_dense + l) =varargout(i,l,j)+1;
        end
    end
    neurons_local{j} = neuron_centers;
    neuron_centers=zeros(0,2);
    for i = 1:num_dense
        for l = 1:num_dense
            this_likelihood = likelihood(i,l,j);
            for ip = max(1, i-x_range) : min(num_dense, i+x_range)
                for lp = max(1, l-y_range) : min(num_dense, l+y_range)
                    dist = (x_dense(i)-x_dense(ip) )^2+ (y_dense(l)-y_dense(lp) )^2;
                        if dist < R
                            if this_likelihood < likelihood(ip,lp,j)
                                this_likelihood = likelihood(ip,lp,j);
                            end
                        end
                end
            end
            if this_likelihood < likelihood_thresholds(threshold_level)
                this_likelihood=0;
            end
            
            if this_likelihood == likelihood(i,l,j)
                neuron_centers = [neuron_centers; [x_dense(i) y_dense(l)] ];
            end
            
        end
    end
    neurons_regional{j} = neuron_centers;
end


