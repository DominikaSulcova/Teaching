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
%           3) resting-state with eyes open / closed
%                   - 1.5 mins each
%                   - at the beginning and in the middle of the session 
% 
% data:     - pre-processed EEG recordings + extracted PSD
% 
% script:   -
%           - 
%           
% output:   1) PAFPain_info 
%           2) PAFPain_data
%           3) PAFPain_measures

%% pararms - ALWAYS RUN AT THE BEGINNING OF THE SESSION
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % MATLAB toolboxes
folder.data = uigetdir(pwd, 'Coose the input folder');          % processed EEG data --> this should be the folder 'PAFPain' at the external harddrive 'PHYSIOLOGIE'
folder.output = uigetdir(pwd, 'Choose the output folder');      % output folder --> local folder with figures, output files, export tables...
cd(folder.output)

% output
study = 'PAFPain';
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

%% extract subject info & pre-processed variables
% ----- section input -----
params.peak = {'N1' 'N2' 'P2'};
% ------------------------- 
% update output structures
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

% load NLEP_info 
if exist('NLEP_info') ~= 1
    load(sprintf('%s\\NLEP_output.mat', folder.data), 'NLEP_info') 
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

%% load and plot PSD
% ----- section input -----
params.log_val = 'off';
params.colours = [0.9216    0.1490    0.1490;
    0.0745    0.6235    1.0000;
    1.0000    0.4784    0.8000;
    0.2588    0.7216    0.0275]; 
% ------------------------- 
% load NLEP_data
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

% save data and continue
PAFPain_data.PSD = data_RSEEG;
save(output_file, 'PAFPain_data', '-append')
clear params a b c d output_vars data visual fig screen_size

%% extract alpha measures per channel 
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

%% compute ICA at alpha frequency
% ----- section input -----
params.prefix = '';
params.bandpass = [8 12];
params.ICA_comp = 7;
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% prepare letswave-friendly dataset
for a = 1:length(data_RSEEG(subject_idx).dataset)
    % prepare header
    dataset(a).header = header;


end

clear a

%% extract alpha measures per component

%% export to excel

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
%             'lims' 'label'
% ========================================================================= 
% define defaults
lims = [-5 5];

% check for varargins
if ~isempty(varargin)
    % limits
    i = find(strcmpi(varargin, 'lims'));
    if ~isempty(i) 
        lims = varargin{i + 1};
    end
end

% plot topography
topoplot(data, chanlocs, 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off', 'maplimits', lims)
set(gca,'color', [1 1 1]);

% add colorbar
colormap(jet)
cb = colorbar;
set(cb, 'Position', [0.15, 0.05, 0.7, 0.03]);
set(cb, 'Orientation', 'horizontal');
end
