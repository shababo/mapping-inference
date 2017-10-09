%% Parameters in the design stage
% Design parameters
n_spots_per_trial = 4;
% Need to run sims to check how these parameters affect the results
n_replicates=1; % conduct two replicates for each trial
K_undefined=5; % each cell appears approximately 10*2 times
K_disconnected=5; % each cell appears approximately 10*2 times
K_connected=10; % each cell appears approximately 10*2 times

single_spot_threshold=9; % switch to single spot stimulation (this can be a function of n_spots_per_trial
trial_max=2000;
% threshold for the group movement:
disconnected_threshold = 0.2;disconnected_confirm_threshold = 0.2;
connected_threshold = 0.5;connected_confirm_threshold = 0.5;
change_threshold=0.01; % for potentially connected cells (since the estimated gamma are important in this case)

% Initialize the five cell groups
undefined_cells= cell(0); undefined_cells{1}=ones(n_cell_this_plane,1);%A
potentially_disconnected_cells= cell(0); potentially_disconnected_cells{1}=zeros(n_cell_this_plane,1);%B
dead_cells= cell(0); dead_cells{1}=zeros(n_cell_this_plane,1);%D
potentially_connected_cells= cell(0); potentially_connected_cells{1}=zeros(n_cell_this_plane,1);%C
alive_cells= cell(0);alive_cells{1}=zeros(n_cell_this_plane,1);%E


% Prior distribution
prior_pi0=0.8;

iter=1;
mpp_undefined=cell(0);trials_locations_undefined=cell(0);trials_powers_undefined=cell(0);
mpp_disconnected=cell(0);trials_locations_disconnected=cell(0);trials_powers_disconnected=cell(0);
mpp_connected=cell(0);trials_locations_connected=cell(0);trials_powers_connected=cell(0);
designs_undefined=[];designs_connected=[];designs_disconnected=[];
outputs_undefined=[];outputs_connected=[];outputs_disconnected=[];


% Initialize tuning parameters in the VI
C_threshold = 0.01;maxit=1000;
S=200;epsilon=0.01;eta_logit=0;
background_rt=background_rate*time_max; % raw probability of firing within a trial

eta_beta=0.05;

% Whether to output plots during the experiment
visualized = 0;

n_trials=0;


gamma_estimates = 0.5*ones(n_cell_this_plane,1);% for drawing samples (not really used)
prob_weight=0;

id_continue=1;% an indicator


% lklh_func=@calculate_likelihood_sum_bernoulli; % likelihood function is
% specificed when fitting the working model

stim_threshold = 10;

% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;

% Used in the random designs
id_notconnected=false;
connected=true;
%  loc_to_cell_nuclei is from get_stim_locations

% Initialize storage
mean_gamma_current=zeros(length(related_cell_list),1);
mean_gain_current=gain_template*ones(length(related_cell_list),1);
gamma_path=zeros(length(related_cell_list),1);
gain_path=zeros(length(related_cell_list),1);
var_gamma_path=zeros(length(related_cell_list),1);

change_threshold = 0.05;
prob_lower_bound = 0.01;


new_rule = false;