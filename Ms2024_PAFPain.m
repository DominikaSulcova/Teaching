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
figure_counter = figure_counter + 1;
saveas(fig, sprintf('%s\\figures\\%s_PSD_original.png', folder.output, PAFPain_info(subject_idx).ID))


clear a b output_vars data fig screen_size

%% extract alpha measures per channel 

%% compute ICA

%% extract alpha measures per component

%% export to excel

%% functions
function plot_PSD(visual, varargin)
% =========================================================================
% Reloads pre-processed EEG data of a single subject for following 
% processing steps. 
% Input:    - visual = a structure with x, y, SD
%           - varargins: 'log_val', 'colours', 'x_lim', 'y_lim', 'shading', 'labels'
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
    if ~isempty(i) strcmp(varargin{i + 1}, 'on')
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