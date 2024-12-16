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
%           - 
%           
% output:   - BetaPain_info     --> MATLAB structure gathering dataset
%                                   information
%           - BetaPain_measures --> MATLAB structure gathering all measures
%                                   intended for statistical analysis
%           - BetaPain_data_SXX --> MATLAB structure gathering all newly
%                                   produced datasets for each subject
%           - excel table with all evaluated variables

%% 1) parameters - ALWAYS RUN AT THE BEGINNING OF THE SESSION
fprintf('section 1:\n')

% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % MATLAB toolboxes
folder.raw = uigetdir(pwd, 'Choose the input folder');              % raw data --> this should be the folder 'BetaPain' at the external harddrive 'PHYSIOLOGIE'
folder.processed = uigetdir(pwd, 'Choose the data folder');         % processed data --> local folder 
folder.output = uigetdir(pwd, 'Choose the output folder');          % output folder --> local folder with figures, output files, export files
cd(folder.processed)

% output
study = 'BetaPain';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

% load the info structure
fprintf('loading the info structure...\n')
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

% current participant
prompt = {'subject number:'};
dlgtitle = 'subject';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
subject_idx = str2num(input{1,1});
clear prompt dlgtitle dims definput input
fprintf('section 1 finished.\n')

% addpath(genpath([folder.toolbox '\letswave 7']));
% letswave

%% 2) participant & session info
fprintf('section 2:\n')

% encode info
fprintf('encoding subject & session info...\n')
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
clear session_info output_vars
fprintf('section 2 finished.\n')

%% 3) load & pre-process EEG data
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};
params.suffix = {'ep' 'ds' 'dc'};
params.eventcode = {'TMS'};
params.epoch = [-3.000 -0.005];
params.downsample = 20;
% ------------------------- 
fprintf('section 3:\n')

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
    dlgtitle = sprintf('subject %d - %s session: dataset', subject_idx, params.condition{a});
    dims = [1 60];
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
    dataset(a).condition = params.condition{a};
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
        if a == 1 && d == 1
            BetaPain_info(subject_idx).EEG.processing(1).process = sprintf('electrode coordinates assigned');
            BetaPain_info(subject_idx).EEG.processing(end).params.layout = sprintf('standard 10-20-cap81');
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % identify eventcode
        fprintf('checking events... ') 
        if d == 1
            params.eventcode_old = unique({lwdata.header.events.code});
            prompt = 'Which eventcode corresponds to the TMS stimulus?';
            definput = '';
            for e = 1:length(params.eventcode_old)
                if e == length(params.eventcode_old)
                    definput = [definput params.eventcode_old{e}];
                else
                    definput = [definput params.eventcode_old{e} '\' ];
                end
            end
            dlgtitle = sprintf('subject %d - %s session: eventcodes', subject_idx, params.condition{a});
            dims = [1 40];
            params.eventcode_old = inputdlg(prompt,dlgtitle,dims,{definput});
            clear prompt dlgtitle dims definput 
        else
            if ~ismember(params.eventcode_old, unique({lwdata.header.events.code}))
                params.eventcode_old1 = params.eventcode_old;
                params.eventcode_old = unique({lwdata.header.events.code});
                prompt = sprintf(['Eventcode ''%s'' was not found in this dataset.\n' ...
                    'Which eventcode corresponds to the TMS stimulus?'], params.eventcode_old1{1});
                definput = '';
                for e = 1:length(params.eventcode_old)
                    if e == length(params.eventcode_old)
                        definput = [definput params.eventcode_old{e}];
                    else
                        definput = [definput params.eventcode_old{e} '\' ];
                    end
                end
                dlgtitle = sprintf('%s', datanames{d});
                dims = [1 40];
                params.eventcode_old = inputdlg(prompt,dlgtitle,dims,{definput});
                clear prompt dlgtitle dims definput 
            end
        end

        % re-label and filter events
        event_idx = logical([]);
        for e = 1:length(lwdata.header.events)
            if strcmp(lwdata.header.events(e).code, params.eventcode_old{1})
                lwdata.header.events(e).code = params.eventcode{1};
                event_idx(e) = false; 
            else
                event_idx(e) = true; 
            end
        end
        lwdata.header.events(event_idx) = [];

        % remove first stimulus if the latency corresponds
        if lwdata.header.events(2).latency - lwdata.header.events(1).latency > 9.8
            lwdata.header.events(1) = [];    
        end
        fprintf('%d TMS stimuli found.\n', length(lwdata.header.events))

        % extract pre-stimulus epochs
        fprintf('epoching ... ')
        option = struct('event_labels', params.eventcode, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
            'x_duration', params.epoch(2)-params.epoch(1), 'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            BetaPain_info(subject_idx).EEG.processing(end+1).process = sprintf('segmented to pre-stimulus epochs');
            BetaPain_info(subject_idx).EEG.processing(end).params.limits = params.epoch;
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{1};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
            epoch_idx = length(BetaPain_info(subject_idx).EEG.processing);
        end
        BetaPain_info(subject_idx).EEG.processing(epoch_idx).params.epochs((a-1)*length(params.timepoint) + d) = size(lwdata.data, 1);

        % downsample 
        fprintf('downsampling... ')
        option = struct('x_dsratio', params.downsample, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_downsample.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            BetaPain_info(subject_idx).EEG.processing(end+1).process = sprintf('downsampled');
            BetaPain_info(subject_idx).EEG.processing(end).params.ratio = params.downsample;
            BetaPain_info(subject_idx).EEG.processing(end).params.fs_orig = 1/lwdata.header.xstep * params.downsample;
            BetaPain_info(subject_idx).EEG.processing(end).params.fs_final = 1/lwdata.header.xstep;
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{2};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % remove DC + linear detrend
        fprintf('removing DC and applying linear detrend.\n')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            BetaPain_info(subject_idx).EEG.processing(end+1).process = sprintf('DC + linear detrend on ERP epochs');
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{3};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).raw(d).header = lwdata.header;
        dataset(a).raw(d).data = lwdata.data; 
    end
    fprintf('Done.\n\n')

    % concatenate timepoints and save for letswave
    fprintf('concatenating blocks and saving... ') 
    for b = 1:length(params.timepoint)
        % identify blocks, check their number
        concat_idx = [];
        for d = 1:length(dataset(a).raw)
            if contains(dataset(a).raw(d).header.name, params.timepoint{b}) 
                concat_idx(end + 1) = d;
            end
        end
        if length(concat_idx) > 2
        elseif length(concat_idx) == 1
            fprintf('ATTENTION: Only 1 block was found for the %s dataset.\n', params.timepoint{b})
        elseif length(concat_idx) == 0
            fprintf('ATTENTION: No %s dataset was found!\n', params.timepoint{b})
            continue
        end

        % merge epochs of all identified datasets
        header = dataset(a).raw(concat_idx(1)).header;
        data = dataset(a).raw(concat_idx(1)).data;
        header.name = header.name(1:end-3); 
        if length(concat_idx) == 2
            for e = 1:size(dataset(a).raw(concat_idx(2)).data, 1)
                data(end+1, :, :, :, :, :) = dataset(a).raw(concat_idx(2)).data(e, :, :, :, :, :);
            end
        end
        header.datasize = size(data);

        % save for letswave
        save([header.name '.lw6'], 'header')
        save([header.name '.mat'], 'data')

        % update dataset
        dataset(a).processed(b).header = header;
        dataset(a).processed(b).data = data;
    end
    fprintf('Done.\n')
    fprintf('\n')
end

% save and continue
save(output_file, 'BetaPain_info','-append')
clear params a b c d e data2import data_idx file_idx filename datanames option lwdata event_idx epoch_idx concat_idx data header
fprintf('section 3 finished.\n')

%% 4) visual inspection
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.prefix = 'dc ds ep';
params.suffix = {'no_mastoid' 'bandpass_broad' 'visual'};
params.bandpass = [1, 45];
% ------------------------- 
fprintf('section 4:\n')

% update output 
load(output_file, 'BetaPain_info')

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% load dataset if needed
if exist('dataset') ~= 1
    data2load = dir(sprintf('%s*%s*', params.prefix, BetaPain_info(subject_idx).ID));
    if length(data2load) == length(params.condition) * length(params.timepoint) * 2
        dataset = reload_dataset(data2load, params.condition, 'processed');
    else
        error(sprintf('ERROR: Wrong number of datasets (%d) found in the directory!', length(data2load)/2))
    end
end

% prepare for viewing
fprintf('exporting for visualization: ')
for a = 1:length(dataset)
    for b = 1:length(dataset(a).processed)
        % select dataset
        lwdata.header = dataset(a).processed(b).header;
        lwdata.data = dataset(a).processed(b).data;

        % remove mastoids
        channel_all = {lwdata.header.chanlocs.labels};
        channel_mask = cellfun(@(x) strcmp(x, 'M1') || strcmp(x, 'M2'), channel_all);
        params.channels2keep = channel_all(~channel_mask);
        option = struct('type', 'channel', 'items', {params.channels2keep}, 'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_selection.get_lwdata(lwdata, option);
        if a == 1 && b == 1
            BetaPain_info(subject_idx).EEG.processing(end + 1).process = 'mastoid channels removed';
            BetaPain_info(subject_idx).EEG.processing(end).params.channels_kept = params.channels2keep;
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{1};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).processed(b).header = lwdata.header;
        dataset(a).processed(b).data = lwdata.data;

        % broad bandpass
        option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2), 'low_cutoff', params.bandpass(1), ... 
            'filter_order', 4, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);

        % save for letswave
        filename = sprintf('%s %s %s %s', params.suffix{3}, BetaPain_info(subject_idx).ID, params.condition{a}, params.timepoint{b});
        header = lwdata.header;
        data = lwdata.data;
        save([filename '.lw6'], 'header')        
        save([filename '.mat'], 'data')
        fprintf('. ')
    end
end
fprintf('done.\n\n')

% save and continue
save(output_file, 'BetaPain_info','-append')
clear params a b data2load lwdata option channel_all channel_mask filename

% open letswave for visual inspection
fprintf('section 4 finished.\nplease perform the visual inspection!\n\n')
letswave

%% 5) compute ICA at beta frequency
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.prefix = 'dc ds ep';
params.suffix = {'no_mastoid' 'reref' 'bandpass_beta' 'ica'};
params.interp_chans = 4;
params.bandpass = [13 30];
params.ICA_comp = 25;
% -------------------------
fprintf('section 5:\n')

% update output 
load(output_file, 'BetaPain_info')

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% load dataset if needed
if exist('dataset') ~= 1
    data2load = dir(sprintf('%s*%s*', params.prefix, BetaPain_info(subject_idx).ID));
    if length(data2load) == length(params.condition) * length(params.timepoint) * 2
        % load the data
        dataset = reload_dataset(data2load, params.condition, 'processed');

        % remove the mastoids
        for a = 1:length(dataset)
            for b = 1:length(dataset(a).processed)
                % select dataset
                lwdata.header = dataset(a).processed(b).header;
                lwdata.data = dataset(a).processed(b).data;

                % remove mastoids
                channel_all = {lwdata.header.chanlocs.labels};
                channel_mask = cellfun(@(x) strcmp(x, 'M1') || strcmp(x, 'M2'), channel_all);
                params.channels2keep = channel_all(~channel_mask);
                option = struct('type', 'channel', 'items', {params.channels2keep}, 'suffix', params.suffix{1}, 'is_save', 0);
                lwdata = FLW_selection.get_lwdata(lwdata, option);

                % update dataset
                dataset(a).processed(b).header = lwdata.header;
                dataset(a).processed(b).data = lwdata.data;
            end 
        end

        % encode if necessary
        encode = true;
        for m = 1:length(BetaPain_info(subject_idx).EEG.processing)
            if contains(BetaPain_info(subject_idx).EEG.processing(m).process, 'mastoid')
                encode = false;
            end
        end
        if encode
            BetaPain_info(subject_idx).EEG.processing(end + 1).process = 'mastoid channels removed';
            BetaPain_info(subject_idx).EEG.processing(end).params.channels_kept = params.channels2keep;
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{1};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end
    else
        error(sprintf('ERROR: Wrong number of datasets (%d) found in the directory!', length(data2load)/2))
    end
end

% interpolate channels if needed
params.labels = {dataset(1).processed(1).header.chanlocs.labels};
prompt = {sprintf('%s session:', params.condition{1}) sprintf('%s session:', params.condition{2})};
dlgtitle = 'channel interpolation';
dims = [1 170];
definput = {strjoin(params.labels,' ') strjoin(params.labels,' ')};
answer = inputdlg(prompt,dlgtitle,dims,definput);
for a = 1:length(answer)
    if ~isempty(answer{a})
        % identify channels to interpolate
        chans2interpolate = split(answer{a}, ' ');

        % interpolate if needed
        for c = 1:length(chans2interpolate)
            if ~isempty(chans2interpolate{c})
                % provide update
                fprintf('%s session: interpolating channel %s\n', params.condition{a}, chans2interpolate{c})

                % indentify the channel to interpolate
                chan_n = find(strcmp(params.labels, chans2interpolate{c}));

                % calculate distances with other electrodes
                chan_dist = -ones(length(dataset(1).processed(1).header.chanlocs), 1);
                for b = setdiff(1:length(dataset(1).processed(1).header.chanlocs), chan_n)
                    if dataset(1).processed(1).header.chanlocs(b).topo_enabled == 1
                        chan_dist(b) = sqrt((dataset(1).processed(1).header.chanlocs(b).X - dataset(1).processed(1).header.chanlocs(chan_n).X)^2 + ...
                            (dataset(1).processed(1).header.chanlocs(b).Y - dataset(1).processed(1).header.chanlocs(chan_n).Y)^2 + ...
                            (dataset(1).processed(1).header.chanlocs(b).Z - dataset(1).processed(1).header.chanlocs(chan_n).Z)^2);
                    end
                end
                chan_dist((chan_dist==-1)) = max(chan_dist);
                [~,chan_dist] = sort(chan_dist);

                % identify neighbouring channels
                chan_dist = chan_dist(1:params.interp_chans);
                chans2use = params.labels;
                chans2use = chans2use(chan_dist);

                % cycle through all datasets
                for d = 1:length(dataset(a).processed)
                    % select data
                    lwdata.header = dataset(a).processed(d).header;
                    lwdata.data = dataset(a).processed(d).data;
        
                    % interpolate using the neighboring electrodes
                    option = struct('channel_to_interpolate', chans2interpolate{c}, 'channels_for_interpolation_list', {chans2use}, ...
                        'suffix', '', 'is_save', 0);
                    lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);
        
                    % update dataset
                    dataset(a).processed(d).header = lwdata.header;
                    dataset(a).processed(d).data = lwdata.data;  
                end
                
                % encode
                if c == 1
                    BetaPain_info(subject_idx).EEG.processing(end+1).process = sprintf('bad channels interpolated');
                    BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
                    encode = length(BetaPain_info(subject_idx).EEG.processing);
                end
                BetaPain_info(subject_idx).EEG.processing(encode).params.bad{c} = chans2interpolate{c};
                BetaPain_info(subject_idx).EEG.processing(encode).params.chans_used{c} = strjoin(chans2use, ' ');  
            end
        end
    end
end
clear prompt dlgtitle dims definput answer

% compute ICA matrix and save for letswave
for a = 1:length(dataset)
    fprintf('\n ======================== %s session ========================\n', params.condition{a})

    % pre-process all datasets
    fprintf('pre-processing: dataset ', params.condition{a})
    for b = 1:length(dataset(a).processed)
        fprintf('%s - ', params.timepoint{b})

        % select dataset
        lwdata.header = dataset(a).processed(b).header;
        lwdata.data = dataset(a).processed(b).data;

        % re-reference to common average
        option = struct('reference_list', {params.labels}, 'apply_list', {params.labels}, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_rereference.get_lwdata(lwdata, option);
        if a == 1 && b == 1
            BetaPain_info(subject_idx).EEG.processing(end+1).process = sprintf('re-referenced to common average');
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{2};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % beta bandpass
        option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2), 'low_cutoff', params.bandpass(1), ... 
            'filter_order', 4, 'suffix', params.suffix{3}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
        if a == 1 && b == 1
            BetaPain_info(subject_idx).EEG.processing(end + 1).process = 'bandpass filtered at beta frequency';
            BetaPain_info(subject_idx).EEG.processing(end).params.method = 'Butterworth';
            BetaPain_info(subject_idx).EEG.processing(end).params.order = 4;
            BetaPain_info(subject_idx).EEG.processing(end).params.limits = params.bandpass;
            BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{3};
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).processed(b).header = lwdata.header;
        dataset(a).processed(b).data = lwdata.data;
    end
    fprintf('done.\n')

    % compute ICA and save  
    fprintf('computing ICA matrix: ')
    lwdataset = dataset(a).processed;
    option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', params.ICA_comp, 'suffix', params.suffix{4}, 'is_save', 1);
    lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    fprintf('done.\n')

    % extract ICA parameters
    matrix(a).condition = params.condition{a};
    matrix(a).mix = lwdataset(1).header.history(end).option.mix_matrix;
    matrix(a).unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
    if a == 1
        params.ICA_chanlocs = lwdataset(1).header.chanlocs;
        for i = 1:size(matrix(a).mix, 2)
            params.ICA_labels{i} = ['IC',num2str(i)];
        end
        params.ICA_fs = 1/lwdataset(1).header.xstep;
    end

    % update dataset
    dataset(a).ICA = lwdataset;

    % unmix data
    for d = 1:length(dataset(a).ICA)
        for e = 1:size(dataset(a).ICA(d).data, 1)
            dataset(a).unmixed(d).header = dataset(a).ICA(d).header;
            dataset(a).unmixed(d).data(e, :, 1, 1, 1, :) = matrix(a).unmix * squeeze(dataset(a).ICA(d).data(e, :, 1, 1, 1, :));        
        end
    end
end
fprintf('\n')

% update info structure
BetaPain_info(subject_idx).EEG.processing(end + 1).process = 'ICA matrix computed';
BetaPain_info(subject_idx).EEG.processing(end).params.method = 'Runica';
BetaPain_info(subject_idx).EEG.processing(end).params.components = params.ICA_comp;
BetaPain_info(subject_idx).EEG.processing(end).params.chanlocs = params.ICA_chanlocs;
BetaPain_info(subject_idx).EEG.processing(end).params.labels = params.ICA_labels;
BetaPain_info(subject_idx).EEG.processing(end).params.fs = params.ICA_fs;
BetaPain_info(subject_idx).EEG.processing(end).params.matrix = matrix;
BetaPain_info(subject_idx).EEG.processing(end).suffix = params.suffix{4};
BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);

% update data structure
for a = 1:length(dataset)
    % condition
    statement = sprintf('BetaPain_data_%s.beta(a).condition = params.condition{a};', BetaPain_info(subject_idx).ID);
    eval(statement)

    % unmixed
    statement = sprintf('BetaPain_data_%s.beta(a).unmixed = dataset(a).unmixed;', BetaPain_info(subject_idx).ID);
    eval(statement)
end
save(output_file, sprintf('BetaPain_data_%s', BetaPain_info(subject_idx).ID), '-append')

% open letswave if not already open
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    letswave
end

% save and continue
save(output_file, 'BetaPain_info', '-append')
clear params a b c d e f i m data2load encode channel_all channel_mask lwdata option chans2interpolate chan_n chan_dist ...
    lwdataset matrix fig_all open data header output_vars
fprintf('section 5 finished.\nplease select components containing sensorimotor beta activity!\n\n')

%% 6) encode selected ICA components
% ----- section input -----
params.condition = {'pain', 'control'};
% ------------------------- 
fprintf('section 6:\n')

% update output 
load(output_file, 'BetaPain_info')

% encode selected ICA components
for a = 1:length(params.condition)
    % ask for the input
    prompt = {'radial components:', 'tangential components:'};
    dlgtitle = sprintf('ICA - %s session', params.condition{a});
    dims = [1 60];
    definput = {'', ''};
    input = inputdlg(prompt,dlgtitle,dims,definput);

    % encode
    BetaPain_info(subject_idx).EEG.processing(end).params.selected(a).condition = params.condition{a};
    BetaPain_info(subject_idx).EEG.processing(end).params.selected(a).radial = str2num(input{1});
    BetaPain_info(subject_idx).EEG.processing(end).params.selected(a).tangential = str2num(input{2});
end

% save and continue
save(output_file, 'BetaPain_info', '-append')
clear params a input prompt dlgtitle dims definput answer 
fprintf('section 6 finished.\n')

%% 7) compute PSD of ICA components
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.foi_limits = [1, 50];
params.method = 'pwelch';
params.comp_type = {'radial' 'tangential'};
params.colours = [0.9216    0.1490    0.1490;
    0.0745    0.6235    1.0000;
    1.0000    0.4784    0.8000;
    0.2588    0.7216    0.0275]; 
% -------------------------   
fprintf('section 7:\n')

% ask for subject number
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% update output 
clear dataset
fprintf('loading the data structure...\n')
load(output_file, 'BetaPain_info', sprintf('BetaPain_data_%s', BetaPain_info(subject_idx).ID))
   
% % add fieldtrip to the top of search path
% addpath(genpath([folder.toolbox '\fieldtrip']));
    
% cycle through sessions
fprintf('computing power spectrum density: \n')
for c = 1:length(params.condition)
    fprintf('%s session: \n', params.condition{c})

    % check for the condition in data
    statement = sprintf('condition = BetaPain_data_%s.beta(c).condition;', BetaPain_info(subject_idx).ID);
    eval(statement)   
    if strcmp(condition, params.condition{c})      
        
        % cycle through timepoints
        for t = 1:length(params.timepoint)
            fprintf('%s - ', params.timepoint{t})
    
            % select the unmixed data and the header
            statement = sprintf('data = BetaPain_data_%s.beta(c).unmixed(t).data;', BetaPain_info(subject_idx).ID);
            eval(statement)        
            statement = sprintf('header = BetaPain_data_%s.beta(c).unmixed(t).header;', BetaPain_info(subject_idx).ID);
            eval(statement) 
             
            % calculate PSD at individual trials
            PSD = struct;
            for e = 1:size(data, 1)
                % % calculate spectra using IRASA
                % % create a FieldTrip data structure
                % cfg = [];
                % cfg.trial = {squeeze(data(e, :, :))};       
                % cfg.time = {0 : header.xstep : (size(data, 3)-1)*header.xstep};      
                % cfg.label = BetaPain_info(subject_idx).EEG.processing(end).params.labels';  
                % data_trial = ft_datatype_raw(cfg);
                % ft_checkdata(data_trial);
                % 
                % % extract original spectra
                % cfg = [];
                % cfg.output = 'pow';
                % cfg.foilim = params.foi_limits;  
                % cfg.pad = 'nextpow2'; 
                % cfg.method = 'irasa';    
                % cfg.output = 'original';
                % PSD(e).original = ft_freqanalysis(cfg, data_trial);
                % 
                % % extract fractal specra
                % cfg.output = 'fractal';
                % PSD(e).fractal = ft_freqanalysis(cfg, data_trial);
                % 
                % % compute oscillatory spectra
                % cfg = [];
                % cfg.parameter = 'powspctrm';
                % cfg.operation     = 'x2-x1';
                % PSD(e).oscillatory = ft_math(cfg, PSD(e).fractal, PSD(e).original);
                
                % calculate original spectrum using pwelch
                for i = 1:size(data, 2)
                    [PSD(e).pwelch.powspctrm(i, :), PSD(e).pwelch.freq] = pwelch(squeeze(data(e, i, :)), [], [], [], 1/header.xstep);
                end
                PSD(e).pwelch.powspctrm(:, PSD(e).pwelch.freq < 1 | PSD(e).pwelch.freq > 50) = [];
                PSD(e).pwelch.freq(PSD(e).pwelch.freq < 1 | PSD(e).pwelch.freq > 50) = [];
            end

            % % create structures
            % psd_irasa(t).params.method = PSD(1).original.cfg.method;
            % psd_irasa(t).params.limits = PSD(1).original.cfg.foilim;
            % psd_irasa(t).params.pad = PSD(1).original.cfg.pad;
            % psd_irasa(t).params.label = PSD(1).original.label';
            % psd_irasa.params.freq = PSD(1).original.freq;
            psd_pwelch(t).params.method = 'pwelch';
            psd_pwelch(t).params.limits = params.foi_limits;
            psd_pwelch(t).params.label = BetaPain_info(subject_idx).EEG.processing(end).params.labels;
            psd_pwelch(t).params.freq = PSD(1).pwelch.freq;
            for e = 1:length(PSD)
                % psd_irasa(t).original(e, :, :) = PSD(e).original.powspctrm;
                % psd_irasa(t).fractal(e, :, :) = PSD(e).fractal.powspctrm;
                % psd_irasa(t).oscillatory(e, :, :) = PSD(e).oscillatory.powspctrm;
                psd_pwelch(t).original(e, :, :) = PSD(e).pwelch.powspctrm;
            end
                
            % % compute average PSD 
            % psd_avg(t).original.mean = squeeze(mean(psd_irasa(t).original, 1));
            % psd_avg(t).original.std = squeeze(std(psd_irasa(t).original, 0, 1));
            % psd_avg(t).fractal.mean = squeeze(mean(psd_irasa(t).fractal, 1));
            % psd_avg(t).fractal.std = squeeze(std(psd_irasa(t).fractal, 0, 1));
            % psd_avg(t).oscillatory.mean = squeeze(mean(psd_irasa(t).oscillatory, 1));
            % psd_avg(t).oscillatory.std = squeeze(std(psd_irasa(t).oscillatory, 0, 1));
            psd_avg(t).pwelch.mean = squeeze(mean(psd_pwelch(t).original, 1));
            psd_avg(t).pwelch.std = squeeze(std(psd_pwelch(t).original, 0, 1));
        end
    else
        error('ERROR: Condition in BetaPain_data does not match condition in section parameters!')
    end

    % % append to the data structure and save
    % statement = sprintf('BetaPain_data_%s.beta(c).psd_irasa = psd_irasa;', BetaPain_info(subject_idx).ID);
    % eval(statement)
    statement = sprintf('BetaPain_data_%s.beta(c).psd_pwelch = psd_pwelch;', BetaPain_info(subject_idx).ID);
    eval(statement)
    statement = sprintf('BetaPain_data_%s.beta(c).psd_avg = psd_avg;', BetaPain_info(subject_idx).ID);
    eval(statement)

    % save data structure 
    save(output_file, sprintf('BetaPain_data_%s', BetaPain_info(subject_idx).ID), '-append')
    fprintf('done.\n')
end
fprintf('\n')

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% select data
statement = sprintf('beta = BetaPain_data_%s.beta;', BetaPain_info(subject_idx).ID);
eval(statement)

% set plotting parameters
switch params.method
    case 'pwelch' 
        visual.x = beta(1).psd_pwelch(1).params.freq';
    case 'irasa'
        visual.x = beta(1).psd(1).params.freq';
end
y_limits = [];

% plot the output figure
fig = figure(figure_counter);
set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
for c = 1:length(params.condition)
    % select component of interest
    coi{1} = BetaPain_info(subject_idx).EEG.processing(end).params.selected(c).radial;
    coi{2} = BetaPain_info(subject_idx).EEG.processing(end).params.selected(c).tangential;
    
    % plot topoplot of the radial component
    if ~isempty(coi{1})
        subplot(4, 15, (c-1)*2*15 + 1)
        topoplot(BetaPain_info(subject_idx).EEG.processing(end).params.matrix(c).mix(:, coi{1}), BetaPain_info(subject_idx).EEG.processing(end).params.chanlocs,...
        'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color', [1 1 1]);
        title(params.comp_type{1})
    else
        subplot(4, 15, (c-1)*2*15 + 1)
        axis([0 1 0 1]); 
        axis off;         
        text(0.5, 0.5, {'radial component', 'not identified'}, 'HorizontalAlignment', 'center', 'FontSize', 10);
    end

    % plot topoplot of the tangential component
    if ~isempty(coi{2})
        subplot(4, 15, (c-1)*2*15 + 15 + 1)
        topoplot(BetaPain_info(subject_idx).EEG.processing(end).params.matrix(c).mix(:, coi{2}), BetaPain_info(subject_idx).EEG.processing(end).params.chanlocs,...
        'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color', [1 1 1]);
        title(params.comp_type{2})
    else
        subplot(4, 15, (c-1)*2*15 + 15 + 1)
        axis([0 1 0 1]); 
        axis off;         
        text(0.5, 0.5, {'tangential component', 'not identified'}, 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    % cycle through timepoints
    for t = 1:length(params.timepoint)
        % plot the psd
        subplot(4, 15, (c-1)*2*15 + [(t-1)*2 + 2, (t-1)*2 + 3, 15 + (t-1)*2 + 2, 15 + (t-1)*2 + 3]) 
        for a = 1:length(coi)
            if ~isempty(coi{a})
                % select data to plot
                switch params.method
                    case 'pwelch'
                        visual.y = beta(c).psd_avg(t).pwelch.mean(coi{a}, :);
                        visual.sd = beta(c).psd_avg(t).pwelch.std(coi{a}, :);
                    case 'irasa'
                        visual.y = beta(c).psd_avg(t).oscillatory.mean(coi{a}, :);
                        visual.sd = beta(c).psd_avg(t).oscillatory.std(coi{a}, :);
                end
    
                % % shade SD
                % S(a) = fill([visual.x fliplr(visual.x)], [visual.y + visual.sd fliplr(visual.y - visual.sd)], ...
                %     params.colours((c-1)*2 + a, :), 'FaceAlpha', 0.2, 'linestyle', 'none');
                % hold on
    
                % plot mean
                P(a) = plot(visual.x, visual.y, 'Color', params.colours((c-1)*2 + a, :), 'LineWidth', 2.2);
                hold on

                % indicate labels
                comp_idx(a) = true;
            else
                % indicate labels
                comp_idx(a) = false;
            end
        end
        
        % add figure title
        if c == 1 && t == 3
            title(sprintf('subject %d: power of beta oscillations of selected components', subject_idx))        
        end

        % set figure limits
        xlim([8, 35])
        y_limits(c, t, :) = get(gca, 'YLim');

        % axes
        box off;
        ax = gca;
        ax.XAxisLocation = 'bottom';
        ax.YAxisLocation = 'left';
        ax.TickDir = 'out'; 
        
        % other parameters
        if c == 2
            xlabel('frequency (Hz)')
        end
        % ylabel(sprintf('PSD (%sV^2/Hz) %s SD)', char(956), char(177)))
        set(gca, 'FontSize', 12)
        set(gca, 'Layer', 'Top')

        % add legend
        if t == 1 
            if ~isempty(coi{1}) && ~isempty(coi{2})
                legend(P, params.comp_type, 'Location', 'northwest', 'Box', 'off');
            end
        end
    end
end

% set scaling across plots
% sgtitle(sprintf('subject %d: power of beta oscillations of selected components', subject_idx))
for c = 1:length(params.condition)
    for t = 1:length(params.timepoint)
        subplot(4, 15, (c-1)*2*15 + [(t-1)*2 + 2, (t-1)*2 + 3, 15 + (t-1)*2 + 2, 15 + (t-1)*2 + 3])         
        ylim([0 max(max(y_limits(1,:)), max(y_limits(2,:)))])            
        
        % % plot frequency limits 
        % hold on
        % L(1) = line([13, 13], [0 max(y_limits, [], 'all')], 'Color', [0.8   0.8    0.8], 'LineWidth', 1.5, 'LineStyle', ':');
        % L(2) = line([30, 30], [0 max(y_limits, [], 'all')], 'Color', [0.8   0.8    0.8], 'LineWidth', 1.5, 'LineStyle', ':');
    end
end

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\%s_beta_%s.png', folder.output, BetaPain_info(subject_idx).ID, params.method))
figure_counter = figure_counter + 1;

% ask for continuation
answer = questdlg('Do you want to continue with next subject?', 'Continue?', 'YES', 'NO', 'YES'); 
if strcmp(answer, 'YES')
    subject_idx = subject_idx + 1;
end
clear params a c e f s t condition data data_trial header PSD psd_irasa psd_pwelch psd_avg cfg coi visual answer y_limits fig S P L 
fprintf('section 7 finished.\n')

%% functions
function dataset = reload_dataset(data2load, conditions, fieldname)
% =========================================================================
% Reloads pre-processed EEG data of a single subject for following 
% processing steps. 
% Input:    - list of datasets to load
%           - cell array with conditions
%           - fieldname
% =========================================================================  
% initiate output
dataset = struct;

% load all data
dpc = length(data2load)/2/length(conditions);
if mod(dpc, 1) == 0
    for c = 1:length(conditions)
        % note condition
        dataset(c).condition = conditions{c}; 

        % subset header and data files
        header_idx = logical([]);
        data_idx = logical([]);
        for d = 1:length(data2load)
            if contains(data2load(d).name, conditions{c}) 
                if contains(data2load(d).name, 'lw6') 
                    header_idx(d) = true;
                    data_idx(d) = false;
                elseif contains(data2load(d).name, 'mat') 
                    header_idx(d) = false;
                    data_idx(d) = true;
                end
            else
                header_idx(d) = false;
                data_idx(d) = false;
            end
        end
        headers = data2load(header_idx);
        datas = data2load(data_idx);

        % load all dataset for this condition
        if length(datas) == length(headers) && length(datas) == dpc
            for d = 1:dpc
                % load header
                load(sprintf('%s\\%s', headers(d).folder, headers(d).name), '-mat')
                statement = sprintf('dataset(c).%s(d).header = header;', fieldname);
                eval(statement) 
    
                % load data
                load(sprintf('%s\\%s', datas(d).folder, datas(d).name))
                statement = sprintf('dataset(c).%s(d).data = data;', fieldname);
                eval(statement) 
            end
        else
            error('ERROR: Wrong number of available datasets to load! Check manually.')
        end
    end
else
    error('ERROR: Wrong number of available datasets to load! Check manually.')
end
end
