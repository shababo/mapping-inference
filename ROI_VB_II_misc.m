%% A map of delta entropy after Stage I
figure(5)
    entropy_locations = weights_mat_full'*delta_H(cell_list);
    entropy_locations(entropy_locations<0) =0; 
    entropies = scatter(Z_dense(:,1),...
        -Z_dense(:,2),...
        abs(entropy_locations)*400+0.01, 'filled','o');
    set(entropies,'MarkerFaceColor','k');
   alpha(entropies,0.7);
    hold on;
    xlim([20,460]);
    ylim([-900,-400]);
    peaksfound = scatter(Z_dense(idx,1),...
        -Z_dense(idx,2),...
        30, 'filled','o');
    set(peaksfound,'MarkerFaceColor','b');
   alpha(peaksfound,0.7);
  
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
        hold off;
     

saveas(5,'../Figures/Entropy_0.jpg')

%% A map of delta entropy after Stage II
figure(6)
    entropy_locations = weights_mat_full'*delta_H(cell_list);
    entropy_locations(entropy_locations<0) =0; 
    entropies = scatter(Z_dense(:,1),...
        -Z_dense(:,2),...
        abs(entropy_locations)*400+0.01, 'filled','o');
    set(entropies,'MarkerFaceColor','k');
   alpha(entropies,0.7);
    hold on;
    xlim([20,460]);
    ylim([-900,-400]);
    peaksfound = scatter(Z_dense(idx,1),...
        -Z_dense(idx,2),...
        30, 'filled','o');
    set(peaksfound,'MarkerFaceColor','b');
   alpha(peaksfound,0.7);
  
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
        hold off;
     
saveas(6,'../Figures/Entropy_100.jpg')

%% Computating time:
delta_time = time_record(2:end) - time_record(1:(ND-1));
figure(12)
plot(delta_time);
hold on;
xlabel('Number of batches');
ylabel('Computing time (s)');
xlim([0,ND]);
ylim([0,2]);
hold off;

saveas(12,'../Data/Computing_time_design.jpg')

%% A map of final estimates 

figure(7)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25);
    set(temp,'MarkerFaceColor','k');
   alpha(temp,0.8);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

    xlim([20,460]);
    ylim([-900,-400]);

for j = 1:size(amp_related_count_trials,2)
    coef = output_post(j).alpha(2:end);
    coef_thres = 0.9;
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), ...
        -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25,'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.4);
    hold on

end

    temp = scatter(postsyn_position(:,1),...
        -postsyn_position(:,2),...
        182,'filled','o');
set(temp,'MarkerFaceColor','g');
   alpha(temp,1);

   
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
saveas(7,'../Figures/Final_fit.jpg')
