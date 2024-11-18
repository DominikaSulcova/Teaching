%% BC2024: The role of sensorimotor beta in pain-motor interactions
% ------------------------------------------------------------------------
% author:   Dominika Sulcova
%           MSH - Medical School Hamburg
% created:  November 2024   
% student:  Antonia Quakernack
% ------------------------------------------------------------------------
% project:  Sensorimotor beta oscillation have been implicated in 
%           functional inhibition of the motor cortex. This thesis project
%           aims to explore whether pain-induced suppression of motor
%           activity is reflected - and can be predicted - by an increase
%           in sensorimotor beta. 
% 
%           Analyzed dataset was acquired in 2020-2022 at the Institute of
%           Neuroscience at UCLouvain, Brussels. 20 healthy subjects were 
%           included, each participating in two experimental sessions:
%           1) pain session     tonic pain was elicited by capsaicin (2%) 
%                               applied in a patch to the hand dorsum
%           2) control session  patch contained only the vehicle (EtOH)
%           In each session, simple biphasic TMS stimuli (120%rMT) were
%           delivered to the M1 in seven successive blocks (80 stim/block):
%           baseline + 6 post-patch blocks. Both EMG and EEG was recorded 
%           but this stidy will only consider MEPs to evaluate pain-induced 
%           inhibition of cortico-motor excitability. 
% 
% data:     - raw EEG recordings 
%               --> 7 blocks, each split in two continuous EEG recordings
%               --> for extraction of pre-stimulus 'resting state' EEG
%           - MEPs
%               --> already pre-processed, saved in 
%           - pain ratings
%               --> already pre-processed, saved in 
%          - subject information
% 
% script:   - subject and session info encoded to the output structure
%           - raw EEG data imported and pre-processed
%               - channels linked to electrode coordinates 
%               - signal segmented to 2.5s chunks pre-stimulus
%               - downsampled to 1000Hz SR 
%               - DC + linear detrend
%               - datasets saved in letswave format 
%           
% output:   1) BetaPain_info 

%% parameters
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
folder.raw = uigetdir(pwd, 'Coose the input folder');           % raw data --> this should be the folder 'BetaPain' at the external harddrive 'PHYSIOLOGIE'
folder.processed = uigetdir(pwd, 'Coose the data folder');      % processed data --> local folder 
folder.output = uigetdir(pwd, 'Choose the output folder');      % output folder --> local folder with figures, output files, export files
cd(folder.processed)

% output
study = 'BetaPain';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

% current participant
prompt = {'subject number:'};
dlgtitle = 'subject';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
subject_idx = str2num(input{1,1});
clear prompt dlgtitle dims definput input

%% participant & session info
% load the info structure
if exist(output_file) == 2
    output_vars = who('-file', sprintf('%s', output_file));
    if ismember('BetaPain_info', output_vars)
        load(output_file, 'BetaPain_info')
    else
        BetaPain_info = struct;
        save(output_file, 'BetaPain_info','-append')
    end
else
    BetaPain_info = struct;
    save(output_file, 'BetaPain_info')
end

% encode info
prompt = {'age:', 'male:', 'handedness:',...
    'session 1 - date:', 'session 1 - condition:', 'session 1 - hemisphere:', 'session 1 - rMT:',...
    'session 2 - date:', 'session 2 - condition:', 'session 2 - hemisphere:', 'session 2 - rMT:'};
dlgtitle = 'Subject session information';
dims = [1 35];
definput = {'18', '1', '1',...
    '00.00.2020', 'pain/control', 'left/right', '',...
    '00.00.2020', 'pain/control', 'left/right', ''};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% update info structure
if subject_idx < 10
    BetaPain_info(subject_idx).ID = sprintf('S0%d', subject_idx);
else
    BetaPain_info(subject_idx).ID = sprintf('S%d', subject_idx);
end
BetaPain_info(subject_idx).age = str2num(session_info{1});
BetaPain_info(subject_idx).male = str2num(session_info{2});
BetaPain_info(subject_idx).handedness = str2num(session_info{3});
BetaPain_info(subject_idx).session(1).date = session_info{4};
BetaPain_info(subject_idx).session(1).condition = session_info{5};
BetaPain_info(subject_idx).stimulation(1).hemisphere = session_info{6};
BetaPain_info(subject_idx).stimulation(1).rMT = str2num(session_info{7});
BetaPain_info(subject_idx).stimulation(1).intensity = BetaPain_info(subject_idx).stimulation(1).rMT * 1.2;
BetaPain_info(subject_idx).session(2).date = session_info{8};
BetaPain_info(subject_idx).session(2).condition = session_info{9};
BetaPain_info(subject_idx).stimulation(2).hemisphere = session_info{10};
BetaPain_info(subject_idx).stimulation(2).rMT = str2num(session_info{11});
BetaPain_info(subject_idx).stimulation(2).intensity = BetaPain_info(subject_idx).stimulation(2).rMT * 1.2;

% save and continue
save(output_file, 'BetaPain_info','-append')
clear session_info

%% load & pre-process EEG data
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};
params.suffix = {'ep' 'ds' 'dc'};
params.eventcode = {'TMS'};
params.epoch = [-2.505 -0.005];
params.downsample = 20;
% ------------------------- 
% cycle through sessions
for a = 1:length(params.condition)
    % add letswave 6 to the top of search path
    addpath(genpath([folder.toolbox '\letswave 6']));

    % provide update
    fprintf('subject %d - %s session:\n', subject_idx, params.condition{a})

    % confirm expected datasets
    datanames = [];
    for b = 1:length(params.timepoint)
        for c = 1:length(params.block)
            datanames{end + 1} = sprintf('%s %s %s %s %s', study, BetaPain_info(subject_idx).ID, params.condition{a}, params.timepoint{b}, params.block{c});
        end
    end
    prompt = {sprintf(['These are expected datasets:\n' ...
    '1 - %s\n2 - %s\n3 - %s\n4 - %s\n5 - %s\n6 - %s\n7 - %s\n' ...
    '8 - %s\n9 - %s\n10 - %s\n11 - %s\n12 - %s\n13 - %s\n14 - %s\n' ...
    'Is any dataset missing (indicate number)?'], ...
    datanames{1},datanames{2},datanames{3},datanames{4},datanames{5},datanames{6},datanames{7}, ...
    datanames{8},datanames{9},datanames{10},datanames{11},datanames{12},datanames{13},datanames{14})};
    dlgtitle = sprintf('subject %d - %s session:\n', subject_idx, params.condition{a});
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    if ~isempty(input{1})
        datanames(str2num(input{1})) = [];
    end
    clear prompt dlgtitle dims definput input

    % identify import folder 
    for b = 1:length(BetaPain_info(subject_idx).session)
        if strcmp(BetaPain_info(subject_idx).session(b).condition, params.condition{a})
            params.folder = sprintf('%s\\%s\\%s', folder.raw, BetaPain_info(subject_idx).ID, ...
                BetaPain_info(subject_idx).session(b).date);
            params.date = BetaPain_info(subject_idx).session(b).date;
        end
    end

    % check avaiable datasets
    data2import = dir(params.folder);
    if ~isempty(data2import)
        data_idx = logical([]);
        for c = 1:length(data2import)
            if isempty(str2num(data2import(c).name))
                data_idx(c) = true;
            else
                data2import(c).block = str2num(data2import(c).name);
                data_idx(c) = false;
            end
        end
        data2import(data_idx) = [];
        [~, file_idx] = sort([data2import.block]);
        data2import = data2import(file_idx);
        fprintf('%d datasets found in the directory.\n', length(data2import))
    else
        error('ERROR: No datasets found in the directory.\n')
    end

    % check number of datasets
    if length(data2import) ~= length(datanames)
        error('ERROR: This does not correspond to the expected number of datasets (%d)!', length(datanames))
    end

    % load datasets
    fprintf('Loading:\n')
    for d = 1:length(data2import)
        % provide update
        fprintf('%s:', datanames{d})
        
        % encode to the info structure
        if a == 1 && d == 1
            BetaPain_info(subject_idx).EEG.dataset(1).session = params.condition{a};
        else
            BetaPain_info(subject_idx).EEG.dataset(end + 1).session = params.condition{a};
        end
        BetaPain_info(subject_idx).EEG.dataset(end).date = params.date;
        BetaPain_info(subject_idx).EEG.dataset(end).block = data2import(d).block;
        BetaPain_info(subject_idx).EEG.dataset(end).name = datanames{d};
    
        % import the dataset
        [dataset(a).raw(d).header, dataset(a).raw(d).data, ~] = RLW_import_MEGA(data2import(d).folder, data2import(d).block);

        % rename in the header
        dataset(a).raw(d).header.name = datanames{d};
    end  
    fprintf('\n')

    % add letswave 7 to the top of search path
    addpath(genpath([folder.toolbox '\letswave 7']));

    % pre-process and save datasets
    fprintf('Pre-processing:\n')
    for d = 1:length(dataset(a).raw)
        % provide update
        fprintf('%s:\n', datanames{d})

        % select data
        lwdata.header = dataset(a).raw(d).header;
        lwdata.data = dataset(a).raw(d).data; 

        % assign electrode coordinates
        fprintf('assigning electrode coordinates... ')
        option = struct('filepath', sprintf('%s\\letswave 7\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', folder.toolbox), ...
            'suffix', '', 'is_save', 0);
        lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            BetaPain_info(subject_idx).EEG.processing(1).process = sprintf('electrode coordinates assigned');
            BetaPain_info(subject_idx).EEG.processing(end).params.layout = sprintf('standard 10-20-cap81');
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % re-label and filter events
        fprintf('checking events... ')
        event_idx = logical([]);
        for b = 1:length(lwdata.header.events)
        end
        lwdata.header.events(event_idx) = [];

    end
    fprintf('Done.\n')
    fprintf('\n')
end

% save and continue
save(output_file, 'BetaPain_info','-append')
clear params a b c d data2import data_idx file_idx filename datanames lwdata event_idx

%% final pre-processing of EEG data
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};
% ------------------------- 