function p = sample_neuron_positions(num_neurons, boundaries)
%create_synthetic_neuron_field(num_neurons, boundaries, pct_excitatory)
%
% Generate num_neurons random synthetic neurons
% 
% num_neurons: number of neurons
% boundaries: size and number of dimensions for neuron space; e.g. [100 100 50]
%       to place neurons in a 100 x 100 x 50 three-dimensional area,
%       centered at [0 0 0]
% pct_excitatory: ratio of excitatory neurons

p = [unifrnd(boundaries(1,1),boundaries(1,2),[num_neurons 1]) ...
     unifrnd(boundaries(2,1),boundaries(2,2),[num_neurons 1]) ...
     unifrnd(boundaries(3,1),boundaries(3,2),[num_neurons 1])];

% ensure neurons aren't overlapping
good_field = 0;
% while ~good_field
%     dist_matrix = squareform(pdist(p));
%     [i,j] = find(dist_matrix < 15 & dist_matrix ~= 0);
%     if isempty(i)
%         good_field = 1;
%     end
% %     for idx = 1:length(i)
% %         % make sure this isn't a duplicate
% %         if i(idx) < j(idx)
% %             p(i(idx),:) = rand(1,length(boundaries)).*boundaries - .5*boundaries;
% %             break
% %         end
% %     end
%     if ~good_field
%         p = [unifrnd(boundaries(1,1),boundaries(1,2),[num_neurons 1]) ...
%              unifrnd(boundaries(1,1),boundaries(1,2),[num_neurons 1]) ...
%              unifrnd(boundaries(1,1),boundaries(1,2),[num_neurons 1])];
%     end
% end

% double check neurons aren't overlapping
dist_matrix = squareform(pdist(p));
if any(dist_matrix < 10 & dist_matrix ~= 0)
    disp('oops')
end
