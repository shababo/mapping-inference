function varargout = mapping_explorer(varargin)
% MAPPING_EXPLORER MATLAB code for mapping_explorer.fig
%      MAPPING_EXPLORER, by itself, creates a new MAPPING_EXPLORER or raises the existing
%      singleton*.
%
%      H = MAPPING_EXPLORER returns the handle to a new MAPPING_EXPLORER or the handle to
%      the existing singleton*.
%
%      MAPPING_EXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPPING_EXPLORER.M with the given input arguments.
%
%      MAPPING_EXPLORER('Property','Value',...) creates a new MAPPING_EXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mapping_explorer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mapping_explorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mapping_explorer

% Last Modified by GUIDE v2.5 15-Dec-2017 16:29:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mapping_explorer_OpeningFcn, ...
                   'gui_OutputFcn',  @mapping_explorer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mapping_explorer is made visible.
function mapping_explorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mapping_explorer (see VARARGIN)

% Choose default command line output for mapping_explorer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mapping_explorer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mapping_explorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_experiment.
function load_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to load_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiopen('matlab')
handles.data.experiment_setup = experiment_setup;
handles.data.neighbourhoods = neighbourhoods;
handles.data.experiment_queries = experiment_queries;

handles.data.neighbourhood_ID = neighbourhoods(1).neighbourhood_ID;
handles.data.batch_ID = 1;

guidata(hObject,handles)

plot_gui(handles)

function plot_gui(handles)

neighbourhood = handles.data.neighbourhoods(...
    [handles.data.neighbourhoods(:,1).neighbourhood_ID] == handles.data.neighbourhood_ID,...
    end);

% experiment id
set(handles.exp_id_text,'string',['Experiment ID: ' handles.data.experiment_setup.exp_id])
set(handles.neighbourhood_ID_text,'string',['Neighbourhood ID: ' num2str(handles.data.neighbourhood_ID)])
set(handles.batch_ID_text,'string',['Batch ID: ' num2str(handles.data.batch_ID)])

primary_neurons = neighbourhood.neurons(~get_group_inds(neighbourhood,'secondary',handles.data.batch_ID));
neuron_locs = get_rowmat_from_structarray(primary_neurons,'location');
all_neurons = get_rowmat_from_structarray(handles.data.experiment_setup.neurons,'location');
neuron_groups = cell(0);
neuron_colors = zeros(size(neuron_locs));
group_colors={'DarkRed', 'DarkGray', 'ForestGreen' 'BlueViolet' 'Black'};
for i_cell = 1:length(primary_neurons)
    neuron_groups{i_cell} = primary_neurons(i_cell).group_ID{handles.data.batch_ID};
    switch primary_neurons(i_cell).group_ID{handles.data.batch_ID}
        case 'disconnected'
            i_group = 1;
        case 'undefined'
            i_group = 2;
        case 'connected'
            i_group = 3;
        case 'alive'
            i_group = 4;
        case 'secondary'
            i_group = 5;
    end
    neuron_colors(i_cell,:) = rgb(group_colors{i_group});
end


% max proj of stack
axes(handles.stack_axes) 
plot_max_proj_w_locs(neighbourhood.stack,neuron_locs,all_neurons,handles.data.experiment_setup);

% plot parameter estimates from posteriors
properties = {'PR_params','gain_params'};
summary_stat = {'lower_quantile','mean','upper_quantile'};
posteriors = grab_values_from_neurons(handles.data.batch_ID,primary_neurons,properties,summary_stat);

axes(handles.gamma_axes)
scatter(neuron_locs(:,2),neuron_locs(:,1),posteriors.PR_params.mean*100,neuron_colors,'filled')
axis ij
ylim([-152 152])
xlim([-152 152])

axes(handles.gain_axes)
scatter(neuron_locs(:,2),neuron_locs(:,1),posteriors.gain_params.mean*2500,neuron_colors,'filled')
axis ij
ylim([-152 152])
xlim([-152 152])

axes(handles.batch_path_axes)
hold off
for i_cell = 1:length(primary_neurons)
    plot([primary_neurons(i_cell).PR_params.mean],'Color',neuron_colors(i_cell,:))
    hold on
end







% --- Executes on button press in neighbourhood_ID_plus.
function neighbourhood_ID_plus_Callback(hObject, eventdata, handles)
% hObject    handle to neighbourhood_ID_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.neighbourhood_ID = handles.data.neighbourhood_ID + 1;
if handles.data.neighbourhood_ID > size(handles.data.neighbourhoods,1)
    handles.data.neighbourhood_ID = 1;
end

guidata(hObject,handles)
plot_gui(handles)

% --- Executes on button press in neighbourhood_ID_minus.
function neighbourhood_ID_minus_Callback(hObject, eventdata, handles)
% hObject    handle to neighbourhood_ID_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.neighbourhood_ID = handles.data.neighbourhood_ID - 1;
if handles.data.neighbourhood_ID < 1
    handles.data.neighbourhood_ID = size(handles.data.neighbourhoods,1);
end

guidata(hObject,handles)
plot_gui(handles)

% --- Executes on button press in batch_ID_minus.
function batch_ID_minus_Callback(hObject, eventdata, handles)
% hObject    handle to batch_ID_minus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.batch_ID = handles.data.batch_ID - 1;
if handles.data.batch_ID < 1
    handles.data.batch_ID = size(handles.data.neighbourhoods,2);
end

guidata(hObject,handles)
plot_gui(handles)

% --- Executes on button press in batch_ID_plus.
function batch_ID_plus_Callback(hObject, eventdata, handles)
% hObject    handle to batch_ID_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.batch_ID = handles.data.batch_ID + 1;
if handles.data.batch_ID > size(handles.data.neighbourhoods,2)
    handles.data.batch_ID = 1;
end

guidata(hObject,handles)
plot_gui(handles)

% --- Executes on button press in refresh_gui.
function refresh_gui_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_gui(handles)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
