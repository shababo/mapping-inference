addpath(genpath('../../mapping-inference'),genpath('../../odessa-beta-beta'));
%% Run get_experiment_setup:

experiment_setup=get_experiment_setup();
experiment_setup.is_exp = 0;
experiment_setup.enable_user_breaks = 0;


%% Generate cells 

    simulation_setup=get_simulation_setup();
    experiment_setup.neurons=generate_neurons(simulation_setup);
%%
neighbourhoods = create_neighbourhoods_caller(experiment_setup);

if experiment_setup.is_exp
    
    handles.data.neighbourhoods = neighbourhoods;
    acq_gui_data.data.neighbourhoods = handles.data.neighbourhoods;
    guidata(acq_gui,acq_gui_data)
    guidata(hObject,handles)
    exp_data = handles.data; save(handles.data.experiment_setup.fullsavefile,'exp_data')
    [acq_gui, acq_gui_data] = get_acq_gui_data;

else
    neighbourhoods = create_neighbourhoods(experiment_setup);   
end

handles = build_first_batch_stim_all_neighborhoods(hObject,handles,acq_gui,acq_gui_data,experiment_setup);

if experiment_setup.is_exp
    % get info on patched cells while first batches prep
    handles = set_cell1_pos(hObject,eventdata,handles,acq_gui,acq_gui_data,experiment_setup);
    [acq_gui, acq_gui_data] = get_acq_gui_data;

    handles = set_cell2_pos(hObject,eventdata,handles,acq_gui,acq_gui_data,experiment_setup);
    [acq_gui, acq_gui_data] = get_acq_gui_data;

    handles = setup_patches(hObject,eventdata,handles,acq_gui,acq_gui_data,experiment_setup);
    [acq_gui, acq_gui_data] = get_acq_gui_data;
    
    % compute bg rate from data

    set(acq_gui_data.test_pulse,'Value',1)
    set(acq_gui_data.loop,'Value',1)
    set(acq_gui_data.tf_on,'Value',get(handles.tf_flag,'Value'));
    set(acq_gui_data.trigger_seq,'Value',1)
    set(acq_gui_data.loop_count,'String',num2str(1))
else
    % simulate bg rate
    experiment_setup.patched_neuron=struct;
    experiment_setup.patched_neuron.background_rate=1e-4;
    experiment_setup.patched_neuron.cell_type=[];
end
num_map_locations = length(neighbourhoods);

not_terminated = 1;
while not_terminated
    for i = 1:num_map_locations
    
        % check for holograms to do
        
%         if experiemnt_setup.experiment_setup.is_exp
%             % move obj
%             set(handles.thenewx,'String',num2str(handles.data.obj_positions(i,1)))
%             set(handles.thenewy,'String',num2str(handles.data.obj_positions(i,2)))
%             set(handles.thenewz,'String',num2str(handles.data.obj_positions(i,3)))
%             [handles,acq_gui,acq_gui_data] = obj_go_to_Callback(handles.obj_go_to,eventdata,handles);
    
    
    
        %------------------------------------------%
        % Run the designed trials
        
        
        if experiment_setup.is_exp
            do_run_trials = 1;
            if experiment_setup.enable_user_breaks
                choice = questdlg('Run the trials?', ...
                    'Run the trials?', ...
                    'Yes','No','Yes');
                % Handle response
                switch choice
                    case 'Yes'
                        do_run_trials = 1;
                        choice = questdlg('Continue user control?',...
                            'Continue user control?', ...
                            'Yes','No','Yes');
                        % Handle response
                        switch choice
                            case 'Yes'
                                experiment_setup.enable_user_breaks = 1;
                            case 'No'
                                experiment_setup.enable_user_breaks = 0;
                        end
                    case 'No'
                        do_run_trials = 0;
                end
            end

            if do_run_trials
                set(handles.tf_flag,'Value',1)
                set(handles.set_seq_trigger,'Value',0)
                set(handles.target_intensity,'String',user_input_powers)

                [handles, acq_gui, acq_gui_data] = build_seq_groups(hObject, eventdata, handles);
            %     handles = guidata(hObject);
            %     acq_gui_data = get_acq_gui_data();
                max_seq_length = str2double(get(handles.max_seq_length,'String'));
                this_seq = acq_gui_data.data.sequence;
                num_runs = ceil(length(this_seq)/max_seq_length);
                handles.data.start_trial = acq_gui_data.data.sweep_counter + 1;

                for run_i = 1:num_runs

                    this_subseq = this_seq((run_i-1)*max_seq_length+1:min(run_i*max_seq_length,length(this_seq)));
                    time_offset = this_subseq(1).start - 1000;
                    for k = 1:length(this_subseq)
                        this_subseq(k).start = this_subseq(k).start - time_offset;
                    end
                    total_duration = (this_subseq(end).start + this_subseq(end).duration)/1000 + 5;

                    set(acq_gui_data.trial_length,'String',num2str(total_duration + 1.0))
                    acq_gui_data = Acq('trial_length_Callback',acq_gui_data.trial_length,eventdata,acq_gui_data);
                    instruction = struct();
                    instruction.type = 32; %SEND SEQ
                    handles.sequence = this_subseq;
                    instruction.sequence = this_subseq;
                    handles.total_duration = total_duration;
                    instruction.waittime = total_duration + 120;
                    disp('sending instruction...')
                    [return_info,success,handles] = do_instruction_slidebook(instruction,handles);
            %         acq_gui_data = get_acq_gui_data();
            %         acq_gui_data.data.stim_key =  return_info.stim_key;
                    acq_gui_data.data.sequence =  this_subseq;
            %         acq_gui = findobj('Tag','acq_gui');
                    guidata(acq_gui,acq_gui_data)

                    set(acq_gui_data.trial_length,'String',num2str(handles.total_duration + 1.0))
                    acq_gui_data = Acq('trial_length_Callback',acq_gui_data.trial_length,eventdata,acq_gui_data);

                    guidata(hObject,handles)
                    exp_data = handles.data; save(handles.data.experiment_setup.fullsavefile,'exp_data')
            %         guidata(acq_gui,acq_gui_data)

                    acq_gui_data = Acq('run_Callback',acq_gui_data.run,eventdata,acq_gui_data);
                    waitfor(acq_gui_data.run,'String','Start')
                    guidata(acq_gui,acq_gui_data)

                end
            end
        else
            % simulate this batch data
            switch experiment_setup.experiment_type
                case 'simulation'
                    % simulate data 
                    i_neighbourhood=i;
                    experiment_query=generate_psc_data(experiment_query,experiment_setup,neighbourhoods(i_neighbourhood));
                case 'reproduction'
                    % read data from files
            end
            
        end
        
        
        % RUN ONLINE MAPPING PIPELINE HERE
        if experiment_setup.is_exp
            % send to analysis computer
        else
            % run online pipeline with data in RAM
            run_online_pipeline(exp_query_filename)
        end
        
        
        
        % Plot the progress
%         fprintf('Number of trials so far: %d; number of cells killed: %d\n',handles.data.design.n_trials{i}, sum(handles.data.design.dead_cells{i}{handles.data.design.iter}+handles.data.design.alive_cells{i}{handles.data.design.iter}))
        
%         do_cont = 0;
%         choice = questdlg('Continue Plane?', ...
%         'Continue Plane?', ...
%         'Yes','No','Yes');
%         % Handle response
%         switch choice
%         case 'Yes'
%             do_cont = 1;
%         case 'No'
%             do_cont = 0;
%         end
% 
%         if ~do_cont   
%             handles.data.design.id_continue{i} = 0;
%         end
    end
end    

exp_data = handles.data; save(handles.data.experiment_setup.fullsavefile,'exp_data')

set(handles.close_socket_check,'Value',1);
instruction.type = 00;
instruction.string = 'done';
[return_info,success,handles] = do_instruction_slidebook(instruction,handles);
guidata(hObject,handles)    