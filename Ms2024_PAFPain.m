%% Ms2024: The relationship between sensorimotor PAF and pain 
% ------------------------------------------------------------------------
% author:   Dominika Sulcova
%           MSH - Medical School Hamburg
% created:  November 2024   
% student:  Lea Matthies
% ------------------------------------------------------------------------
% project:  The peak latency of alpha oscillation have been previously
%           linked to individual pain sensitivity. However, there are many
%           sources of alpha activity in the brain and it is likely that
%           only some of these are related to pain processing. This study 
%           aims to investigate, whether sensorimotor alpha can predict 
%           laser-evoked brain responses (LEPs) and/or the magnitude of 
%           elicited pain. We hypothesize that a relationship between PAF 
%           and pain should be specific to sensorimotor areas (as compared 
%           to visual areas) and to pain-related brain response (as 
%           compared to non-painful somatosensory-evoked potentials, SEPs).     
% 
%           Analyzed dataset was acquired in 2024 in the pain research lab
%           at the Medical School Hamburg (MSH). 45 healthy subjects were 
%           included, each participating in a single experimental session.
%           During the experiment, 63-cahnnel EEG was recorded during:         
%           1) painful laser stimulation AND innocuous electric stimulation
%                   - areas: both hands / both feet / a hand and a foot
%                   - 2 blocks of 30 stimuli per each area
%           2) resting-state with eyes open / closed
%                   - 1.5 mins each
%                   - at the beginning and in the middle of the session 
% 
% data:     - pre-processed EEG recordings: RS-EEG, LEPs, SEPs 
%           - full-spectrum PSD extracted from all RS-EEG datasets
%           - pain ratings provided by subjects
% 
% script:   - extracts subject info & pre-processed variables
%           - loads and plots full-spectrum PSD
%           - extracts alpha measures per channel, plots topo distributions
%           - computes ICA at alpha frequency and plots all components
%           
% output:   - PAFPain_info      --> MATLAB structure gathering dataset
%                                   information
%           - PAFPain_data      --> MATLAB structure gathering all newly
%                                   produced datasets
%           - PAFPain_measures  --> MATLAB structure gathering all measures
%                                   intended for statistical analysis
%           - excel table with all evaluated variables

%% 1) pararms - ALWAYS RUN AT THE BEGINNING OF THE SESSION
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
folder.data = uigetdir(pwd, 'Coose the input folder');          % processed EEG data 
folder.output = uigetdir(pwd, 'Choose the output folder');      % output folder --> local folder with figures, output files, export tables...
cd(folder.output)

% output
study = 'PAFPain';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

% load output structures
fprintf('loading output structures...\n')
if exist(output_file) == 2 
    % check available variables
    output_vars = who('-file', sprintf('%s', output_file));

    % info structure
    if ismember('PAFPain_info', output_vars)
        load(output_file, 'PAFPain_info')
    else
        PAFPain_info = struct;
        save(output_file, 'PAFPain_info','-append')
    end

    % data structure
    if ismember('PAFPain_data', output_vars)
        load(output_file, 'PAFPain_data')
    else
        PAFPain_info = struct;
        save(output_file, 'PAFPain_data','-append')
    end

    % measures structure
    if ismember('PAFPain_measures', output_vars)
        load(output_file, 'PAFPain_measures')
    else
        PAFPain_info = struct;
        save(output_file, 'PAFPain_measures','-append')
    end
else
    PAFPain_info = struct;
    PAFPain_data = struct;
    PAFPain_measures = struct;
    save(output_file, 'PAFPain_info', 'PAFPain_data', 'PAFPain_measures')
end

% current participant
prompt = {'subject number:'};
dlgtitle = 'subject';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
subject_idx = str2num(input{1,1});
clear prompt dlgtitle dims definput input output_vars
fprintf('section finished.\n')

%% 2) extract subject info & pre-processed variables
% ----- section input -----
params.peak = {'N1' 'N2' 'P2'};
% ------------------------- 
% load NLEP_info 
if exist('NLEP_info') ~= 1
    load('NLEP_output.mat', 'NLEP_info') 
end

% extract subject & session information
fprintf('\nsubject %d:\n', subject_idx)
fprintf('extracting subject information...\n')
PAFPain_info(subject_idx).ID = NLEP_info.single_subject(subject_idx).ID;
PAFPain_info(subject_idx).age = NLEP_info.single_subject(subject_idx).age;
PAFPain_info(subject_idx).male = NLEP_info.single_subject(subject_idx).male;
PAFPain_info(subject_idx).handedness = NLEP_info.single_subject(subject_idx).handedness;
PAFPain_info(subject_idx).area{1} = NLEP_info.single_subject(subject_idx).condition{1}(1:4);
PAFPain_info(subject_idx).area{2} = NLEP_info.single_subject(subject_idx).condition{2}(1:4);
PAFPain_info(subject_idx).side{1} = NLEP_info.single_subject(subject_idx).condition{1}(6:end);
PAFPain_info(subject_idx).side{2} = NLEP_info.single_subject(subject_idx).condition{2}(6:end);

% load NLEP_measures 
if exist('NLEP_measures') ~= 1
    load(sprintf('%s\\NLEP_output.mat', folder.data), 'NLEP_measures') 
end

% extract LEP measures
fprintf('extracting LEP measures...\n')
PAFPain_measures(subject_idx).LEP.conditions = NLEP_measures(subject_idx).LEP_avg.conditions; 
PAFPain_measures(subject_idx).LEP.peaks = params.peak; 
% average LEP values
for c = 1:length(PAFPain_measures(subject_idx).LEP.conditions)
    for p = 1:length(PAFPain_measures(subject_idx).LEP.peaks)
        if strcmp(PAFPain_measures(subject_idx).LEP.peaks{p}, 'N1')
            % N1: import average measures extracted from the ICA-filtered data
            statement = sprintf('PAFPain_measures(subject_idx).LEP.average.amplitude(c, p) = mean(NLEP_measures(subject_idx).LEP_avg.%s.ICA_filtered.amplitude(c, :));', params.peak{p});
            eval(statement)
            statement = sprintf('PAFPain_measures(subject_idx).LEP.average.latency(c, p) = mean(NLEP_measures(subject_idx).LEP_avg.%s.ICA_filtered.latency(c, :));', params.peak{p});
            eval(statement)
        else
            % N2, P2: import average measures extracted from raw data
            statement = sprintf('PAFPain_measures(subject_idx).LEP.average.amplitude(c, p) = mean(NLEP_measures(subject_idx).LEP_avg.%s.raw.amplitude(c, :));', params.peak{p});
            eval(statement)
            statement = sprintf('PAFPain_measures(subject_idx).LEP.average.latency(c, p) = mean(NLEP_measures(subject_idx).LEP_avg.%s.raw.latency(c, :));', params.peak{p});
            eval(statement)
        end
    end
end
% single-trial LEP values
for c = 1:length(PAFPain_measures(subject_idx).LEP.conditions)
    % check that conditions match
    if strcmp(NLEP_measures(subject_idx).LEP_ST.conditions{c}, PAFPain_measures(subject_idx).LEP.conditions{c})
        for p = 1:length(PAFPain_measures(subject_idx).LEP.peaks)
            PAFPain_measures(subject_idx).LEP.single_trial.amplitude(c, p, :) = squeeze(NLEP_measures(subject_idx).LEP_ST.amplitude(p, c, :));
            PAFPain_measures(subject_idx).LEP.single_trial.latency(c, p, :) = squeeze(NLEP_measures(subject_idx).LEP_ST.latency(p, c, :));
        end
    else
        error('ERROR: conditions do not match!')
    end
end

% extract pain ratings
fprintf('extracting pain ratings...\n')
PAFPain_measures(subject_idx).pain.conditions = PAFPain_measures(subject_idx).LEP.conditions;
for c = 1:length(PAFPain_measures(subject_idx).pain.conditions)    
    PAFPain_measures(subject_idx).pain.ratings{c} = [];
    for d = find(contains(NLEP_measures(subject_idx).conditions, PAFPain_measures(subject_idx).pain.conditions{c}))'
        PAFPain_measures(subject_idx).pain.ratings{c}(end + 1: end + length(NLEP_measures(subject_idx).pain(d, :))) = NLEP_measures(subject_idx).pain(d, :);
    end
end
fprintf('done.\n')

% save and continue
save(output_file, 'PAFPain_info', 'PAFPain_measures', '-append')
clear params c d p output_vars statement
fprintf('section finished.\n')

%% 3) load and plot PSD
% ----- section input -----
params.log_val = 'off';
params.colours = [0.9216    0.1490    0.1490;
    0.0745    0.6235    1.0000;
    1.0000    0.4784    0.8000;
    0.2588    0.7216    0.0275]; 
% ------------------------- 
% load NLEP_data
if exist('data_RSEEG') ~= 1
    output_vars = who('-file', 'NLEP_output.mat');
    for a = 1:length(output_vars)
        if contains(output_vars{a}, 'data') 
            load('NLEP_output.mat', output_vars{a})
            if evalin('base', sprintf('isstruct(%s)', output_vars{a}))
            else
                evalin('base', sprintf('clear %s', output_vars{a}))
            end
        end
    end
    data_RSEEG = NLEP_data.RSEEG;
    data_RSEEG(1:35) = NLEP_data_1to35.RSEEG(1:35);
    data_RSEEG = rmfield(data_RSEEG, "PSD_avg");
    clear NLEP_data NLEP_data_1to35
end

% extract mean values for plotting
visual.x = data_RSEEG(subject_idx).freq;
data = [];
for b = 1:length(data_RSEEG(subject_idx).PSD_st)
    data(1, end + 1 : end + size(data_RSEEG(subject_idx).PSD_st(b).original, 1), :) = squeeze(mean(data_RSEEG(subject_idx).PSD_st(b).original, 2));
    data(2, end + 1 : end + size(data_RSEEG(subject_idx).PSD_st(b).fractal, 1), :) = squeeze(mean(data_RSEEG(subject_idx).PSD_st(b).fractal, 2));
end
visual.y{1} = mean(squeeze(data(1, :, :)), 1);
visual.y{2} = mean(squeeze(data(2, :, :)), 1);
visual.SD{1} = std(squeeze(data(1, :, :)), 0, 1);
visual.SD{2} = std(squeeze(data(2, :, :)), 0, 1);

% plot mean PSD
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3) / 2.5, screen_size(4) / 2])
plot_PSD(visual,'log_val', 'on', 'x_lim', [5, 80], ...
    'colours', params.colours(1:length(visual.y), :), 'labels', {'original spectrum',  'aperiodic component'})

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\%s_PSD_original.png', folder.output, PAFPain_info(subject_idx).ID))
figure_counter = figure_counter + 1;

% get average signals across conditions
data_conditions{1} = []; data_conditions{2} = []; data_conditions{3} = []; data_conditions{4} = [];
for c = 1:length(data_RSEEG(subject_idx).dataset)
    if contains(data_RSEEG(subject_idx).dataset{c}, 'open')
        data_conditions{1}(end + 1 : end + size(data_RSEEG(subject_idx).PSD_st(c).oscillatory, 1), :, :) = data_RSEEG(subject_idx).PSD_st(c).oscillatory;
    elseif contains(data_RSEEG(subject_idx).dataset{c}, 'closed')
        data_conditions{2}(end + 1 : end + size(data_RSEEG(subject_idx).PSD_st(c).oscillatory, 1), :, :) = data_RSEEG(subject_idx).PSD_st(c).oscillatory;
    elseif contains(data_RSEEG(subject_idx).dataset{c}, sprintf('%s %s', PAFPain_info(subject_idx).area{1}, PAFPain_info(subject_idx).side{1}))
        data_conditions{3}(end + 1 : end + size(data_RSEEG(subject_idx).PSD_st(c).oscillatory, 1), :, :) = data_RSEEG(subject_idx).PSD_st(c).oscillatory;
    elseif contains(data_RSEEG(subject_idx).dataset{c}, sprintf('%s %s', PAFPain_info(subject_idx).area{2}, PAFPain_info(subject_idx).side{2}))
        data_conditions{4}(end + 1 : end + size(data_RSEEG(subject_idx).PSD_st(c).oscillatory, 1), :, :) = data_RSEEG(subject_idx).PSD_st(c).oscillatory;
    end
end
for d = 1:length(data_conditions)
    visual.y{d} = squeeze(mean(data_conditions{d}, [1,2]))';
    visual.SD{d} = squeeze(std(squeeze(mean(data_conditions{d}, 2)), 0, 1));
end

% plot oscillatory PSD
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3) / 2.5, screen_size(4) / 2])
plot_PSD(visual,'log_val', params.log_val, 'x_lim', [5, 45], 'shading', 'off', ...
    'colours', params.colours(1:length(visual.y), :), ...
    'labels', {'RS-EEG eyes open', 'RS-EEG eyes closed', ...
    sprintf('%s %s', PAFPain_info(subject_idx).area{1}, PAFPain_info(subject_idx).side{1}), ...
    sprintf('%s %s', PAFPain_info(subject_idx).area{2}, PAFPain_info(subject_idx).side{2})})

% save figure and update counter
saveas(fig, sprintf('%s\\figures\\%s_PSD_oscillatory.png', folder.output, PAFPain_info(subject_idx).ID))
figure_counter = figure_counter + 1;

% identify conditions
params.conditions = {'open' 'closed' ...
    sprintf('%s %s', PAFPain_info(subject_idx).area{1}, PAFPain_info(subject_idx).side{1})...
    sprintf('%s %s', PAFPain_info(subject_idx).area{2}, PAFPain_info(subject_idx).side{2})};

% save average data to data structure
PAFPain_data(subject_idx).PSD_avg = struct([]);
for c = 1:length(params.conditions)
    % identify and select dataset(s)
    data_idx = contains(data_RSEEG(subject_idx).dataset, params.conditions{c}, 'IgnoreCase', true); 
    data_codnition.dataset = data_RSEEG(subject_idx).dataset(data_idx);
    data_codnition.psd = data_RSEEG(subject_idx).PSD_st(data_idx);

    % cycle through datsets
    for d = 1:length(data_codnition.dataset)
        % encode condition
        PAFPain_data(subject_idx).PSD_avg(end+1).condition = params.conditions{c};

        % encode dataset name
        PAFPain_data(subject_idx).PSD_avg(end).dataset = data_codnition.dataset{d};

        % encode frequencies
        PAFPain_data(subject_idx).PSD_avg(end).freqs = data_RSEEG(subject_idx).freq;

        % save data averaged per codition
        PAFPain_data(subject_idx).PSD_avg(end).original = squeeze(mean(data_codnition.psd(d).original, 1));
        PAFPain_data(subject_idx).PSD_avg(end).fractal = squeeze(mean(data_codnition.psd(d).fractal, 1));
        PAFPain_data(subject_idx).PSD_avg(end).oscillatory = squeeze(mean(data_codnition.psd(d).oscillatory, 1));
    end
end
save(output_file, 'PAFPain_data', '-append')
clear params a b c d output_vars data visual fig screen_size data_idx data_codnition
fprintf('section finished.\n')

%% 4) extract alpha measures per channel 
% ----- section input -----
params.FOI = [7  15];
params.measures = {'PAF' 'CAF' 'amplitude'};
% -------------------------
% determine frequencies of interest
freq_idx = data_RSEEG(subject_idx).freq >= params.FOI(1) & data_RSEEG(subject_idx).freq <= params.FOI(2);
freq_psd = data_RSEEG(subject_idx).freq(freq_idx);

% cycle through all datasets/trials/electrodes
for a = 1:length(data_RSEEG(subject_idx).dataset)
    for b = 1:size(data_RSEEG(subject_idx).PSD_st(a).oscillatory, 1)
        for c = 1:size(data_RSEEG(subject_idx).PSD_st(a).oscillatory, 2)
            % select data at frequencies of interest            
            psd_trial = squeeze(data_RSEEG(subject_idx).PSD_st(a).oscillatory(b, c, freq_idx))';

            % determine PAF and CAF
            PAF{a}(b, c) = freq_psd(find(psd_trial == max(psd_trial)));

            % determine centroid frequency
            if min(psd_trial) < 0
                CAF{a}(b, c) = sum(freq_psd .* (psd_trial - min(psd_trial) + 1)) / sum(psd_trial - min(psd_trial) + 1);
            else
                CAF{a}(b, c) = sum(freq_psd .* psd_trial) / sum(psd_trial);
            end

            % determine amplitude
            amplitude{a}(b, c) = max(psd_trial);            
        end
    end
end

% split data into conditions
labels_conditions = {'open' 'closed' ...
    sprintf('%s %s', PAFPain_info(subject_idx).area{1}, PAFPain_info(subject_idx).side{1}) ...
    sprintf('%s %s', PAFPain_info(subject_idx).area{2}, PAFPain_info(subject_idx).side{2})};
PAF_conditions{1} = []; PAF_conditions{2} = []; PAF_conditions{3} = []; PAF_conditions{4} = [];
CAF_conditions{1} = []; CAF_conditions{2} = []; CAF_conditions{3} = []; CAF_conditions{4} = [];
amplitude_conditions{1} = []; amplitude_conditions{2} = []; amplitude_conditions{3} = []; amplitude_conditions{4} = [];
for c = 1:length(data_RSEEG(subject_idx).dataset)
    if contains(data_RSEEG(subject_idx).dataset{c}, labels_conditions{1})
        PAF_conditions{1}(end + 1 : end + size(PAF{c}, 1), :) = PAF{c};
        CAF_conditions{1}(end + 1 : end + size(CAF{c}, 1), :) = CAF{c};
        amplitude_conditions{1}(end + 1 : end + size(amplitude{c}, 1), :) = amplitude{c};
    elseif contains(data_RSEEG(subject_idx).dataset{c}, labels_conditions{2})
        PAF_conditions{2}(end + 1 : end + size(PAF{c}, 1), :) = PAF{c};
        CAF_conditions{2}(end + 1 : end + size(CAF{c}, 1), :) = CAF{c};
        amplitude_conditions{2}(end + 1 : end + size(amplitude{c}, 1), :) = amplitude{c};
    elseif contains(data_RSEEG(subject_idx).dataset{c}, labels_conditions{3})
        PAF_conditions{3}(end + 1 : end + size(PAF{c}, 1), :) = PAF{c};
        CAF_conditions{3}(end + 1 : end + size(CAF{c}, 1), :) = CAF{c};
        amplitude_conditions{3}(end + 1 : end + size(amplitude{c}, 1), :) = amplitude{c};
    elseif contains(data_RSEEG(subject_idx).dataset{c}, labels_conditions{4})
        PAF_conditions{4}(end + 1 : end + size(PAF{c}, 1), :) = PAF{c};
        CAF_conditions{4}(end + 1 : end + size(CAF{c}, 1), :) = CAF{c};
        amplitude_conditions{4}(end + 1 : end + size(amplitude{c}, 1), :) = amplitude{c};
    end
end

% save to measures structure
PAFPain_measures(subject_idx).alpha.conditions = labels_conditions;
PAFPain_measures(subject_idx).alpha.PAF = PAF_conditions;
PAFPain_measures(subject_idx).alpha.CAF = CAF_conditions;
PAFPain_measures(subject_idx).alpha.amplitude = amplitude_conditions;

% load default header
load('dataset_default.lw6', '-mat')

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% plot topoplots with mean values
fig_CAF = figure('Position', [100, 100, 700, 550]);
fig_amplitude = figure('Position', [150, 150, 700, 550]);
for d = 1:length(labels_conditions)
    % extract mean measures
    CAF_plot = mean(CAF_conditions{d}, 1);
    amplitude_plot = mean(amplitude_conditions{d}, 1);

    % plot CAF topoplot
    figure(fig_CAF)
    subplot(2, 2, d)
    plot_topo(CAF_plot, header.chanlocs, 'lims', [8 12])
    title(labels_conditions{d})
    if d == 1
        sgtitle('CAF distribution')
    end
    hold on

    % plot CAF topoplot
    figure(fig_amplitude)
    subplot(2, 2, d)
    plot_topo(amplitude_plot, header.chanlocs, 'lims', [0 10])
    if d == 1
        sgtitle('alpha amplitude distribution')
    end
    title(labels_conditions{d})
    hold on
end

% save figures
saveas(fig_CAF, sprintf('%s\\figures\\%s_CAF_topo.png', folder.output, PAFPain_info(subject_idx).ID))
saveas(fig_amplitude, sprintf('%s\\figures\\%s_amplitude_topo.png', folder.output, PAFPain_info(subject_idx).ID))

% save output structures and continue
save(output_file, 'PAFPain_measures', '-append')
clear params a b c d freq_idx freq_psd psd_trial ...
    PAF PAF_conditions CAF CAF_plot CAF_conditions amplitude amplitude_plot amplitude_conditions ...
    fig_amplitude fig_CAF labels_conditions 
fprintf('section finished.\n')

%% 5) compute ICA at alpha frequency
% ----- section input -----
params.prefix = {'icfilt ica_all chunked' 'icfilt ica_all RS'};
params.suffix = {'alpha' 'ica'};
params.bandpass = [7 13];
params.ICA_comp = 8;
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% identify conditions
params.conditions = {'open' 'closed' ...
    sprintf('%s %s', PAFPain_info(subject_idx).area{1}, PAFPain_info(subject_idx).side{1})...
    sprintf('%s %s', PAFPain_info(subject_idx).area{2}, PAFPain_info(subject_idx).side{2})};

% look for available RS-EEG data
data2load = dir(sprintf('%s\\*%s\\%s*%s*', folder.data, PAFPain_info(subject_idx).ID, params.prefix{1}, PAFPain_info(subject_idx).ID));
data2load = [data2load; dir(sprintf('%s\\*%s\\%s*%s*', folder.data, PAFPain_info(subject_idx).ID, params.prefix{2}, PAFPain_info(subject_idx).ID))];

% load data into a letswave-friendly dataset
dataset = struct([]);
dataset_counter = 1;
if length(data2load) == length(data_RSEEG(subject_idx).dataset)*2
    % update
    fprintf('loading the datasets: ')

    % cycle through conditions
    for a = 1:length(params.conditions)
        % identify the datasets
        dataset2load = data2load(contains({data2load.name}, params.conditions{a}, 'IgnoreCase', true)); 

        % cycle through the datsets
        for b = 1:length(dataset2load)                
            if contains(dataset2load(b).name, 'lw6')
                % update
                fprintf('%d ... ', dataset_counter)
                dataset_counter = dataset_counter + 1;

                % append condition
                dataset(end+1).condition = params.conditions{a};

                % append name
                dataset(end).name = dataset2load(b).name;

                % load the header
                load(sprintf('%s\\%s', dataset2load(b).folder, dataset2load(b).name), '-mat');
                dataset(end).header = header;
            elseif contains(dataset2load(b).name, 'mat')
                % load the data
                load(sprintf('%s\\%s', dataset2load(b).folder, dataset2load(b).name));
                dataset(end).data = data;
            end
        end
    end
else
    error('ERROR: Wrong number of datasets (%d) found in the directory!', length(data2load))
end
fprintf('done.\n\n')

% filter in the alpha band
fprintf('applying alpha bandpass filter:\ndataset ')
for c = 1:length(dataset)
    % update
    fprintf('%d ...', c)

    % select dataset
    lwdata.header = dataset(c).header;
    lwdata.data = dataset(c).data;

    % alpha bandpass
    option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2), 'low_cutoff', params.bandpass(1), ... 
        'filter_order', 4, 'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
    if c == 1
        PAFPain_info(subject_idx).ICA(1).process = 'bandpass filtered at alpha frequency';
        PAFPain_info(subject_idx).ICA(1).params.method = 'Butterworth';
        PAFPain_info(subject_idx).ICA(1).params.order = 4;
        PAFPain_info(subject_idx).ICA(1).params.limits = params.bandpass;
        PAFPain_info(subject_idx).ICA(1).suffix = params.suffix{1};
        PAFPain_info(subject_idx).ICA(1).date = sprintf('%s', date);
    end

    % update dataset
    dataset(c).header = lwdata.header;
    dataset(c).data = lwdata.data;
end
fprintf('done.\n\n')

% subset lwdataset
fileds2remove = {'condition' 'name'};
lwdataset = rmfield(dataset, fileds2remove);

% compute ICA and save  
fprintf('computing ICA matrix: \n')
option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', params.ICA_comp, 'suffix', params.suffix{2}, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
fprintf('done.\n\n')

% update dataset
for d = 1:length(dataset)
    dataset(d).header = lwdataset(d).header;
    dataset(d).data = lwdataset(d).data;
end

% extract ICA parameters
matrix.mix = dataset(1).header.history(end).option.mix_matrix;
matrix.unmix = dataset(1).header.history(end).option.unmix_matrix;   
for i = 1:size(matrix.mix, 2)
    params.ICA_labels{i} = ['IC', num2str(i)];
end
params.ICA_fs = 1/dataset(1).header.xstep;

% encode ICA parameters
PAFPain_info(subject_idx).ICA(2).process = 'ICA matrix computed';
PAFPain_info(subject_idx).ICA(2).params.method = 'Runica';
PAFPain_info(subject_idx).ICA(2).params.components = params.ICA_comp;
PAFPain_info(subject_idx).ICA(2).params.chanlocs = dataset(1).header.chanlocs;
PAFPain_info(subject_idx).ICA(2).params.labels = params.ICA_labels;
PAFPain_info(subject_idx).ICA(2).params.fs = params.ICA_fs;
PAFPain_info(subject_idx).ICA(2).params.matrix_mix = matrix.mix;
PAFPain_info(subject_idx).ICA(2).params.matrix_unmix = matrix.unmix;
PAFPain_info(subject_idx).ICA(2).suffix = params.suffix{2};
PAFPain_info(subject_idx).ICA(2).date = sprintf('%s', date);

% unmix ICA components and save to data structure
for d = 1:length(dataset)
    % encode condition
    PAFPain_data(subject_idx).ICA(d).condition = dataset(d).condition;

    % save the header
    PAFPain_data(subject_idx).ICA(d).header = dataset(d).header;
    
    % save unmixed data
    for e = 1:size(dataset(d).data, 1)        
        PAFPain_data(subject_idx).ICA(d).data_unmixed(e, :, 1, 1, 1, :) = matrix.unmix * squeeze(dataset(d).data(e, :, 1, 1, 1, :));        
    end
end
save(output_file, 'PAFPain_data', '-append')

% calculate MRS of unmixed data
PAFPain_measures(subject_idx).ICA.conditions = params.conditions;
data_conditions{1} = []; data_conditions{2} = []; data_conditions{3} = []; data_conditions{4} = [];
for d = 1:length(PAFPain_data(subject_idx).ICA)
    % extract MRS at single-trial level
    data = [];
    for e = 1:size(PAFPain_data(subject_idx).ICA(d).data_unmixed, 1)
        for i = 1:size(PAFPain_data(subject_idx).ICA(d).data_unmixed, 2)
            data(e, i) = sqrt(mean(squeeze(PAFPain_data(subject_idx).ICA(d).data_unmixed(e, i, 1, 1, 1, :)).^2));
        end
    end

    % append data
    if contains(PAFPain_data(subject_idx).ICA(d).condition, PAFPain_measures(subject_idx).ICA.conditions{1})
        data_conditions{1}(end + 1 : end + size(data, 1), :) = data;
    elseif contains(PAFPain_data(subject_idx).ICA(d).condition, PAFPain_measures(subject_idx).ICA.conditions{2})
        data_conditions{2}(end + 1 : end + size(data, 1), :) = data;
    elseif contains(PAFPain_data(subject_idx).ICA(d).condition, PAFPain_measures(subject_idx).ICA.conditions{3})
        data_conditions{3}(end + 1 : end + size(data, 1), :) = data;
    elseif contains(PAFPain_data(subject_idx).ICA(d).condition, PAFPain_measures(subject_idx).ICA.conditions{4})
        data_conditions{4}(end + 1 : end + size(data, 1), :) = data;
    end
end
PAFPain_measures(subject_idx).ICA.MRS = data_conditions;

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% prepare a structure for plotting
visual.components = PAFPain_info(subject_idx).ICA(2).params.labels;  
visual.labels_visual = PAFPain_measures(subject_idx).ICA.conditions([1,2]);  
visual.labels_sensory = PAFPain_measures(subject_idx).ICA.conditions([3,4]); 
visual.chanlocs = PAFPain_info(subject_idx).ICA(2).params.chanlocs;  

% plot IC topographies and MRS differences across eyes open/closed datasets
fig_visual = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig_visual, 'Position', [screen_size(3)/5, 1, 3*screen_size(3)/5, 9*screen_size(4)/10])
for i = 1:length(visual.components)
    % plot IC topography 
    subplot(ceil(length(visual.components)/2), 4, (i-1)*2 + 1);
    visual.topo = double(PAFPain_info(subject_idx).ICA(2).params.matrix_mix (:, i)');
    plot_topo(visual.topo, visual.chanlocs, 'lims', [-3, 3], 'colorbar', 'off')
    title(visual.components{i})

    % calculate mean MRS
    visual.MRS(1) = mean(data_conditions{1}(:, i));
    visual.MRS(2) = mean(data_conditions{2}(:, i));
    visual.SEM(1) = std(data_conditions{1}(:, i)) / sqrt(size(data_conditions{1}(:, i), 1));
    visual.SEM(2) = std(data_conditions{2}(:, i)) / sqrt(size(data_conditions{2}(:, i), 1));

    % plot barplot
    subplot(ceil(length(visual.components)/2), 4, (i-1)*2 + 2);
    plot_barplot(visual.MRS, visual.SEM, visual.labels_visual, 'y_label', sprintf('mean MRS %s SEM', char(177)))
end
sgtitle('ICs: sensitivity to visual input')

% save figure and update figure counter
saveas(fig_visual, sprintf('%s\\figures\\%s_ICs_visual.png', folder.output, PAFPain_info(subject_idx).ID))
figure_counter = figure_counter + 1;

% plot IC topographies and MRS differences across stimulated areas
fig_sensory = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig_sensory, 'Position', [screen_size(3)/5, 1, 3*screen_size(3)/5, 9*screen_size(4)/10])
for i = 1:length(visual.components)
    % plot IC topography 
    subplot(ceil(length(visual.components)/2), 4, (i-1)*2 + 1);
    visual.topo = double(PAFPain_info(subject_idx).ICA(2).params.matrix_mix (:, i)');
    plot_topo(visual.topo, visual.chanlocs, 'lims', [-3, 3], 'colorbar', 'off')
    title(visual.components{i})

    % calculate mean MRS
    visual.MRS(1) = mean(data_conditions{3}(:, i));
    visual.MRS(2) = mean(data_conditions{4}(:, i));
    visual.SEM(1) = std(data_conditions{3}(:, i)) / sqrt(size(data_conditions{3}(:, i), 1));
    visual.SEM(2) = std(data_conditions{4}(:, i)) / sqrt(size(data_conditions{4}(:, i), 1));

    % plot barplot
    subplot(ceil(length(visual.components)/2), 4, (i-1)*2 + 2);
    plot_barplot(visual.MRS, visual.SEM, visual.labels_sensory, 'y_label', sprintf('mean MRS %s SEM', char(177)), 'palette', 'winter')
end
sgtitle('ICs: sensitivity to sensory input')

% save figure and update figure counter
saveas(fig_sensory, sprintf('%s\\figures\\%s_ICs_sensory.png', folder.output, PAFPain_info(subject_idx).ID))
figure_counter = figure_counter + 1;

% save the measures structure 
save(output_file, 'PAFPain_measures', '-append')

% categorize components
prompt = {'visual ICs:' 'bilateral sensorimotor ICs:' sprintf('%s ICs:', params.conditions{3}) sprintf('%s ICs:', params.conditions{4})};
dlgtitle = 'categorize ICA components';
dims = [1 50];
definput = {'' '' '' ''};
input = inputdlg(prompt,dlgtitle,dims,definput);

replace(params.conditions{3}, ' ', '_')

% encode to the info structure and save
PAFPain_info(subject_idx).ICA(3).process = 'ICs categorized';
PAFPain_info(subject_idx).ICA(3).params.method = 'visual inspection';
PAFPain_info(subject_idx).ICA(3).params.visual = str2num(input{1});
PAFPain_info(subject_idx).ICA(3).params.sm_bilateral = str2num(input{2});
statement = sprintf('PAFPain_info(subject_idx).ICA(3).params.sm_%s = str2num(input{3});', replace(params.conditions{3}, ' ', '_'));
eval(statement)
statement = sprintf('PAFPain_info(subject_idx).ICA(3).params.sm_%s = str2num(input{4});', replace(params.conditions{4}, ' ', '_'));
eval(statement)
PAFPain_info(subject_idx).ICA(3).date = sprintf('%s', date);
save(output_file, 'PAFPain_info', '-append')

% ask for continuation
answer = questdlg('Do you want to continue with next subject?', 'Continue?', 'YES', 'NO', 'YES'); 
if strcmp(answer, 'YES')
    subject_idx = subject_idx + 1;
end
clear params a b c d e i data2load dataset2load data header dataset_counter ...
    lwdata option lwdataset fileds2remove matrix data_conditions ...
    visual fig_visual fig_sensory screen_size prompt dlgtitle dims definput input
fprintf('section finished.\n')

%% 6) extract alpha measures per component
% ----- section input -----
params.NFFT = 1000;
params.FOI_save = [5, 20];
params.FOI = [7  13];
% -------------------------
% cycle though all subjects 
for subject_idx = 2:length(PAFPain_data)
    % provide update
    fprintf('subject %d:\n', subject_idx)

    % extract single-trial alpha measures
    fprintf('processing single-trials... ')
    for a = 1:length(PAFPain_measures(subject_idx).ICA.conditions)
        % launch the output 
        PAF{a} = []; CAF{a} = []; amplitude{a} = []; 

        % cycle through datasets
        for b = 1:length(PAFPain_data(subject_idx).ICA)
            % encode component labels
            for d = 1:size(PAFPain_data(subject_idx).ICA(b).data_unmixed, 2)
                PAFPain_data(subject_idx).ICA(b).components{d} = sprintf('IC%d', d); 
            end
                    
            % if it belongs to the condition
            if strcmp(PAFPain_measures(subject_idx).ICA.conditions{a}, PAFPain_data(subject_idx).ICA(b).condition)
                for c = 1:size(PAFPain_data(subject_idx).ICA(b).data_unmixed, 1)
                    % cycle through components
                    for d = 1:size(PAFPain_data(subject_idx).ICA(b).data_unmixed, 2)
                        % calculate PSD
                        [PSD.powspctrm, PSD.freq] = pwelch(squeeze(PAFPain_data(subject_idx).ICA(b).data_unmixed(c, d, 1, 1, 1, :))', ...
                            [], [], params.NFFT, 1/PAFPain_data(subject_idx).ICA(b).header.xstep);

                        % save for frequencies of interest
                        if c == 1 && d == 1
                            PAFPain_data(subject_idx).ICA(b).freq = PSD.freq(PSD.freq > params.FOI_save(1) & PSD.freq < params.FOI_save(2)); 
                        end
                        PAFPain_data(subject_idx).ICA(b).PSD(c, d, :) = PSD.powspctrm(PSD.freq > params.FOI_save(1) & PSD.freq < params.FOI_save(2));

                        % crop around alpha band
                        freq_idx = PSD.freq >= params.FOI(1) & PSD.freq <= params.FOI(2);
                        PSD.freq = PSD.freq(freq_idx);
                        PSD.powspctrm = PSD.powspctrm(freq_idx);
                        
                        % determine PAF
                        if d == 1
                            PAF{a}(end + 1, d) = PSD.freq(find(PSD.powspctrm == max(PSD.powspctrm)));
                        else
                            PAF{a}(end, d) = PSD.freq(find(PSD.powspctrm == max(PSD.powspctrm)));
                        end

                        % determine CAF
                        if min(PSD.powspctrm) < 0
                            if d == 1
                                CAF{a}(end + 1, d) = sum(PSD.freq .* (PSD.powspctrm - min(PSD.powspctrm) + 1)) / sum(PSD.powspctrm - min(PSD.powspctrm) + 1);
                            else
                                CAF{a}(end, d) = sum(PSD.freq .* (PSD.powspctrm - min(PSD.powspctrm) + 1)) / sum(PSD.powspctrm - min(PSD.powspctrm) + 1);
                            end
                        else
                            if d == 1
                                CAF{a}(end + 1, d) = sum(PSD.freq .* PSD.powspctrm) / sum(PSD.powspctrm);
                            else
                                CAF{a}(end, d) = sum(PSD.freq .* PSD.powspctrm) / sum(PSD.powspctrm);
                            end
                        end

                        % determine amplitude
                        if d == 1
                            amplitude{a}(end + 1, d) = max(PSD.powspctrm);  
                        else
                            amplitude{a}(end, d) = max(PSD.powspctrm);
                        end
                    end
                end
            end
        end
    end

    % save to the output structure
    PAFPain_measures(subject_idx).ICA.components = PAFPain_data(subject_idx).ICA(b).components; 
    PAFPain_measures(subject_idx).ICA.PAF = PAF;
    PAFPain_measures(subject_idx).ICA.CAF = CAF;
    PAFPain_measures(subject_idx).ICA.amplitude = amplitude;

    % calculate average PSD
    fprintf('processing average PSD... ')
    for b = 1:length(PAFPain_data(subject_idx).ICA)
        PAFPain_data(subject_idx).ICA(b).PSD_avg = squeeze(mean(PAFPain_data(subject_idx).ICA(b).PSD, 1));
    end

    % extract average alpha measures per condition
    PAFPain_measures(subject_idx).ICA_avg.conditions = PAFPain_measures(subject_idx).ICA.conditions;
    PAFPain_measures(subject_idx).ICA_avg.components = PAFPain_data(subject_idx).ICA(b).components; 
    for a = 1:length(PAFPain_measures(subject_idx).ICA.conditions)
        % launch data array
        data{a} = [];

        % get average PSD per condition
        for b = 1:length(PAFPain_data(subject_idx).ICA)
            % if it belongs to the condition
            if strcmp(PAFPain_measures(subject_idx).ICA.conditions{a}, PAFPain_data(subject_idx).ICA(b).condition) 
                % append to data
                data{a}(end + 1 : end + length(PAFPain_data(subject_idx).ICA(b).PSD), :, :) = PAFPain_data(subject_idx).ICA(b).PSD;
            end
        end
        freq_idx = PAFPain_data(subject_idx).ICA(1).freq >= params.FOI(1) & PAFPain_data(subject_idx).ICA(1).freq <= params.FOI(2);
        freq_avg = PAFPain_data(subject_idx).ICA(1).freq(freq_idx)';
        data_avg{a} = squeeze(mean(data{a}, 1)); 
        data_avg{a} = data_avg{a}(:, freq_idx);

        % determine alpha measures
        for d = 1:size(data_avg{a}, 1)
            % PAF
            PAFPain_measures(subject_idx).ICA_avg.PAF{a}(d) = freq_avg(data_avg{a}(d, :) == max(data_avg{a}(d, :)));

            % CAF
            if min(data_avg{a}(d, :)) < 0
                PAFPain_measures(subject_idx).ICA_avg.CAF{a}(d) = sum(freq_avg .* (data_avg{a}(d, :) - min(data_avg{a}(d, :)) + 1)) / sum(data_avg{a}(d, :) - min(data_avg{a}(d, :)) + 1);
            else
                PAFPain_measures(subject_idx).ICA_avg.CAF{a}(d) = sum(freq_avg .* data_avg{a}(d, :)) / sum(data_avg{a}(d, :));
            end

            % amplitude
            PAFPain_measures(subject_idx).ICA_avg.amplitude{a}(d) = max(data_avg{a}(d, :));
        end
    end
    fprintf('\n')
end
fprintf('\ndone.\n')    

% save and continue
save(output_file, 'PAFPain_data', 'PAFPain_measures', '-append') 
clear params a b c d PSD PAF CAF amplitude data freq_idx freq_avg data_avg 
fprintf('section finished.\n')

%% calculate GFP of LEPs and SEPs
%% export to excel
%% plot group values
%% functions
function plot_PSD(visual, varargin)
% =========================================================================
% Plots the PSD.  
% Input:    - visual = a structure with x, y, SD
%           - varargins: 
%             'log_val', 'colours', 'x_lim', 'y_lim', 'shading', 'labels'
% ========================================================================= 
% set defaults
log_val = false;
colours = parula(length(visual.y));  
x_lim = [visual.x(1), visual.x(end)];
shading = true;
alpha = 0.2;
legend_loc = 'northeast';
labels = [];
for a = 1:length(visual.y)
    labels{a} = sprintf('signal %d', a);
end

% check for varargins
if ~isempty(varargin)
    % log input values
    i = find(strcmpi(varargin, 'log_val'));
    if ~isempty(i) && strcmp(varargin{i + 1}, 'on')
        log_val = true;
    end

    % colours
    i = find(strcmpi(varargin, 'colours'));
    if ~isempty(i)
        colours = varargin{i + 1};
    end

    % x limits
    i = find(strcmpi(varargin, 'x_lim'));
    if ~isempty(i)
        x_lim = varargin{i + 1};
    end

    % y limits
    i = find(strcmpi(varargin, 'y_lim'));
    if ~isempty(i)
        y_lim = varargin{i + 1};
    end

    % shading 
    i = find(strcmpi(varargin, 'shading'));
    if ~isempty(i) && strcmp(varargin{i + 1}, 'off')
        shading = false;
    end

    % labels
    i = find(strcmpi(varargin, 'labels'));
    if ~isempty(i)
        labels = varargin{i + 1};       
    end     
end

% prepare data and labels
if log_val
    % log the input values
    visual.x = log10(visual.x);
    for a = 1:length(visual.y)
        visual.y{a} = log10(visual.y{a});
        visual.SD{a} = log10(visual.SD{a});
    end

    % adjust limits
    x_lim = log10(x_lim);
    if exist("y_lim")
    end
    
    % determine legend location
    legend_loc = 'southwest';

    % prepare the axis labels 
    ax_labels.x = 'frequency log10(Hz)'; 
    ax_labels.y = 'PSD log10(μV²/Hz)'; 
else
    % determine legend location
    legend_loc = 'northeast';

    % prepare the axis labels 
    ax_labels.x = 'frequency (Hz)'; 
    ax_labels.y = 'PSD (μV²/Hz)'; 
end

% shade standard deviation
if shading
    for a = 1:length(visual.y) 
        fill([visual.x fliplr(visual.x)], [visual.y{a} + visual.SD{a} fliplr(visual.y{a} - visual.SD{a})], ...
            colours(a, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        hold on
    end
end

% plot mean values
for a = 1:length(visual.y) 
    P(a) = plot(visual.x, visual.y{a}, 'Color', colours(a, :), 'LineWidth', 2.5);
    hold on
end

% deal with limis
xlim(x_lim)
if ~exist('y_lim')
    y_lim = ylim;
end
ylim(y_lim)

% legend
legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
legend('boxoff');

% set other parameters
xlabel(ax_labels.x)
ylabel(ax_labels.y)
set(gca, 'FontSize', 14)
set(gca, 'Layer', 'Top')
end
function plot_topo(data, chanlocs, varargin)
% =========================================================================
% Plots a topoplot.  
% Input:    - data = a vector of measured values 1 x n of channels
%           - chanlocs = a structure with channel locations
%           - varargins: 
%             'lims' 'colorbar'
% ========================================================================= 
% define defaults
lims = [-5 5];
plot_cb = true;

% check for varargins
if ~isempty(varargin)
    % limits
    i = find(strcmpi(varargin, 'lims'));
    if ~isempty(i) 
        lims = varargin{i + 1};
    end

    % colorbar
    i = find(strcmpi(varargin, 'colorbar'));
    if ~isempty(i) && strcmp(varargin{i + 1}, 'off')
        plot_cb = false;
    end
end

% plot topography
topoplot(data, chanlocs, 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off', 'maplimits', lims)
set(gca,'color', [1 1 1]);

% add colorbar
colormap(jet)
if plot_cb
    cb = colorbar;
    set(cb, 'Position', [0.15, 0.05, 0.7, 0.03]);
    set(cb, 'Orientation', 'horizontal');
end
end
function plot_barplot(means, whiskers, labels, varargin)
% =========================================================================
% Plots a simple barplot.  
% Input:    - means = a vector of mean values
%           - whiskers = a vector of SEM or SD
%           - labels = a cell array with group labels
%           - varagins: 'y_label', 'colours', 'palette'
% ========================================================================= 
% set defaults
y_label = sprintf('mean value %s SEM', char(177));
colours = prism(length(means)); 

% check for varargins
if ~isempty(varargin)
    % y axis label
    i = find(strcmpi(varargin, 'y_label'));
    if ~isempty(i) 
        y_label = varargin{i + 1};
    end

    % colours
    i = find(strcmpi(varargin, 'colours'));
    if ~isempty(i) 
        colours = varargin{i + 1};
    end

    % palette
    i = find(strcmpi(varargin, 'palette'));
    if ~isempty(i) 
        statement = sprintf('colours = %s(length(means)); ', varargin{i + 1});
        eval(statement)
    end
end

% plot bar plot
B = bar(means, 'FaceColor', 'flat');  
hold on;

% add one-sided error bars
E = errorbar(1:2, means, zeros(size(whiskers)), whiskers, 'k', 'linestyle', 'none', 'LineWidth', 1.5);

% other parameters
set(gca, 'XTickLabel', labels); 
ylabel(y_label);
B.CData = colours;
end
