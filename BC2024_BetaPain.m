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
% data:     1) raw EEG recordings 
%               --> 7 blocks, each split in two continuous EEG recordings
%               --> for extraction of pre-stimulus 'resting state' EEG
%           2) MEPs
%               --> already pre-processed, saved in 
%           3) pain ratings
%               --> already pre-processed, saved in 
%           4) subject information
% 
% script:
% output:

%% parameters
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
folder.raw = uigetdir(pwd, 'Coose the input folder');           % raw data --> at MSH, this should be the study folder at the V drive
folder.output = uigetdir(pwd, 'Choose the OneDrive folder');    % output folder --> One Drive: figures, loutput file, exports 
cd(folder.output)

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

%% load EEG data
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};
% ------------------------- 
% define datanames
datanames = [];
for a = 1:length(params.condition)
    for b = 1:length(params.timepoint)
        for c = 1:length(params.block)
            datanames{end + 1} = sprintf('%s %s %s %s %s', study, BetaPain_info(subject_idx).ID, params.condition{a}, params.timepoint{b}, params.block{c});
        end
    end
end

% cycle through sessions
for a = 1:length(params.condition)
    % provide update
    fprintf('subject %d - %s session:\n', subject_idx, params.condition{a})

    % confirm expected datasets
    datanames_a = datanames((a-1)*length(params.timepoint) + [1:14]);
    prompt = {sprintf(['These are expected datasets:\n' ...
    '1 - %s\n2 - %s\n3 - %s\n4 - %s\n5 - %s\n6 - %s\n7 - %s\n' ...
    '8 - %s\n9 - %s\n10 - %s\n11 - %s\n12 - %s\n13 - %s\n14 - %s\n' ...
    'Is any dataset missing (indicate number)?'], ...
    datanames_a{1},datanames_a{2},datanames_a{3},datanames_a{4},datanames_a{5},datanames_a{6},datanames_a{7}, ...
    datanames_a{8},datanames_a{9},datanames_a{10},datanames_a{11},datanames_a{12},datanames_a{13},datanames_a{14})};
    dlgtitle = sprintf();
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    clear prompt dlgtitle dims definput input

    % identify import folder 
    for b = 1:length(BetaPain_info(subject_idx).session)
        if strcmp(BetaPain_info(subject_idx).session(b).condition, params.condition{a})
            params.folder = sprintf('%s\\%s\\%s', folder.raw, BetaPain_info(subject_idx).ID, ...
                BetaPain_info(subject_idx).session(b).date);
        end
    end

    % check avaiable datasets
    data2import = dir(params.folder);
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
    fprintf('%d datasets found in the directory. Loading:\n', length(data2import))

    % load datasets
    for d = 1:length(data2import)
        % identify names
        filename = sprintf('%s\\%s', data2import(d).folder, data2import(d).name);
        dataname = datanames{};

        % provide update
        fprintf('%s ...\n', dataname)
        
        % encode the filename to metadata
        block = regexp(file2import(c).name, 'b(\d+)', 'tokens');
        block = str2num(block{1}{1});
        RFSxLASER_info(subject_idx).dataset(block - 1).block = block;
        RFSxLASER_info(subject_idx).dataset(block - 1).name = file2import(c).name;
    
        % import the dataset
        [dataset(subject_idx).raw(block - 1).header, dataset(subject_idx).raw(block - 1).data, ~] = RLW_import_VHDR(filename);
    
        % rename in the header
        dataset(subject_idx).raw(block - 1).header.name = dataname;
    end
   
end

clear params a b c d data2import data_idx file_idx filename dataname