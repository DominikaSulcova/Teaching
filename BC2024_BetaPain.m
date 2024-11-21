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
% output:   1) BetaPain_info 
%           2) BetaPain_measures
%           3) BetaPain_data

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
clear session_info output_vars

%% load & pre-process EEG data
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.block = {'b1', 'b2'};
params.suffix = {'ep' 'ds' 'dc'};
params.eventcode = {'TMS'};
params.epoch = [-3.000 -0.005];
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
    dlgtitle = sprintf('subject %d - %s session: dataset', subject_idx, params.condition{a});
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
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % remove DC + linear detrend
        fprintf('removing DC and applying linear detrend.\n')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if a ==1 && d == 1
            BetaPain_info(subject_idx).EEG.processing(end+1).process = sprintf('DC + linear detrend on ERP epochs');
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).raw(d).header = lwdata.header;
        dataset(a).raw(d).data = lwdata.data; 
    end
    fprintf('Done.\n')
    fprintf('\n')

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

%% visual inspection
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.prefix = 'dc ds ep';
params.suffix = {'no_mastoid' 'bandpass_broad' 'visual'};
params.bandpass = [1, 45];
% ------------------------- 
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
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).processed(b).header = lwdata.header;
        dataset(a).processed(b).data = lwdata.data;

        % broad bandpass
        option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2), 'low_cutoff', params.bandpass(1), ... 
            'filter_order', 4, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
        if a == 1 && b == 1
            BetaPain_info(subject_idx).EEG.processing(end + 1).process = 'bandpass filtered for visualization';
            BetaPain_info(subject_idx).EEG.processing(end).params.method = 'Butterworth';
            BetaPain_info(subject_idx).EEG.processing(end).params.order = 4;
            BetaPain_info(subject_idx).EEG.processing(end).params.limits = params.bandpass;
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

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
letswave

%% extraction of sensorimotor beta oscillations
% ----- section input -----
params.condition = {'pain', 'control'};
params.timepoint = {'baseline', 't1', 't2', 't3', 't4', 't5', 't6'};
params.prefix = 'dc ds ep';
params.suffix = {'no_mastoid' 'bandpass_beta'};
params.interp_chans = 4;
% ------------------------- 
% re-reference to common average
% filter at beta frequency 15 - 30Hz?
% separately for each condition:    - ICA 25 channels
%                                   - select non-artifactual components
%                                     over stimulated hemisphere
%                                   - encode selected components
% extract beta amplitude and PBL:   - frequency decomposition of retained
%                                     signal at each electrode/trial                                  
%                                   - average across electrodes/trials
%                                   - remove aperiodic component
%                                   - calculate peak amplitude and latency                                 

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

% pre-process
fprintf(': ')
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
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end

        % update dataset
        dataset(a).processed(b).header = lwdata.header;
        dataset(a).processed(b).data = lwdata.data;

        % broad bandpass
        option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2), 'low_cutoff', params.bandpass(1), ... 
            'filter_order', 4, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
        if a == 1 && b == 1
            BetaPain_info(subject_idx).EEG.processing(end + 1).process = 'bandpass filtered for visualization';
            BetaPain_info(subject_idx).EEG.processing(end).params.method = 'Butterworth';
            BetaPain_info(subject_idx).EEG.processing(end).params.order = 4;
            BetaPain_info(subject_idx).EEG.processing(end).params.limits = params.bandpass;
            BetaPain_info(subject_idx).EEG.processing(end).date = sprintf('%s', date);
        end
    end
end

% save and continue
save(output_file, 'BetaPain_info','-append')
clear params a b c m data2load encode channel_all channel_mask lwdata option chans2interpolate chan_n chan_dist

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
