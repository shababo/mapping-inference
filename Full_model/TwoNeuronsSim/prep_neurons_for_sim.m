function neurons=prep_neurons_for_sim(neurons,loc_target,z_dictionary,z_shape_function,z_variance)


for i=1:2
   rel_loc = loc_target-neurons(i).location;
   true_shape_func=z_dictionary{neurons(i).shape_index};
   neurons(i).z=struct;

   neurons(i).z.true_received=max(0,true_shape_func(rel_loc));
   neurons(i).z.received=max(0,z_shape_function(rel_loc));
   neurons(i).z.variance=max(0,z_variance(rel_loc));
   
end

