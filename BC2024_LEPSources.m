%% BC2024: Cortical sources of laser-evoked potentials
% ------------------------------------------------------------------------
% author:   Dominika Sulcova
%           MSH - Medical School Hamburg
% created:  November 2024   
% student:  Pauline Kraus
% ------------------------------------------------------------------------
% project:  
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
% data:     
% 
% script:   
%           
% output:   
% 

%% params
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');            % local folder --> MATLAB toolboxes
folder.LEP = uigetdir(pwd, 'Choose the folder with LEP data');          % MSH D-drive folder --> processed LEPs 
folder.output = uigetdir(pwd, 'Choose the output folder');              % output folder --> local folder with figures, output files, export files
cd(folder.output)

% output
study = 'LEPSources';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

%% export selected datasets for BESA
% ----- section input -----
params.subjects = [1	2	4	7	9	10	11	14	15	16	17	18	19	22	23	24	25	26	27	28	29	30	31	32	33	34	35	43	44	45];
params.side = {'right' 'right' 'right' 'left' 'right' 'left' 'left' 'right' 'left' 'left' ...
    'right' 'left' 'left' 'right' 'left' 'left' 'left' 'left' 'right' 'left' ...
    'left' 'right' 'right' 'left' 'right' 'right' 'right' 'left' 'left' 'right'};
params.stimulus = 'LEP';
params.area = 'hand';
params.block = {'b1' 'b2'};
params.prefix = 'bl icfilt ica_all dc ep reref ds notch bandpass dc'; 
% ------------------------- 
% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% load single-trial data
dataset = struct;
fprintf('loading:\nsubject')
for s = 1:length(params.subjects)
    % provide update
    fprintf(' %d -', s)

    % determine subject ID
    if params.subjects(s) < 10
        subject_ID = sprintf('S00%d', params.subjects(s));
    else
        subject_ID = sprintf('S0%d', params.subjects(s));
    end

    % load both blocks
    dataset_block = struct;
    for b = 1:length(params.block)
        % check available files
        data2load = dir(sprintf('%s\\NLEP_%s\\%s %s %s %s %s %s*', ...
            folder.LEP, subject_ID, params.prefix, subject_ID, params.stimulus, params.area, params.side{s}, params.block{b}));
        if length(data2load) == 2
            for d = 1:length(data2load)
                if contains(data2load(d).name, 'lw6')
                    % load header
                    load(sprintf('%s\\%s', data2load(d).folder, data2load(d).name), '-mat')
                    dataset_block(b).header = header;
                elseif contains(data2load(d).name, 'mat')
                    % load data
                    load(sprintf('%s\\%s', data2load(d).folder, data2load(d).name))
                    dataset_block(b).data = data;
                end
            end
        else
            error('ERROR: Incorrect number of files (%d) found in the directory!', length(data2load))
        end
    end

    % concatenate blocks
    header = dataset_block(1).header;
    header.name = header.name(1:end-3);
    data = dataset_block(1).data;
    for e = 1:size(dataset_block(2).data, 1)
        data(end+1, :, :, :, :, :) = dataset_block(2).data(e, :, :, :, :, 1:size(data, 6));
    end
    header.datasize = size(data);
    dataset.original(s).header = header;
    dataset.original(s).data = data;   
end
clear b d e s subject_ID dataset_block data2load data header  
fprintf('\ndone.\n\n')

% prepare flip dictionary
labels = {dataset.original(1).header.chanlocs.labels};
labels_flipped = labels;
for i = 1:length(labels)
    electrode_n = str2num(labels{i}(end));
    if isempty(electrode_n)
    else
        if mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
            a = find(strcmpi(labels, label_new));
            if isempty(a)
            else
                labels_flipped{i} = label_new;
            end
        else                                    % even number --> right hemisphere 
            label_new = labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
            a = find(strcmpi(labels, label_new));
            if isempty(a)
            else
                labels_flipped{i} = label_new;
            end
        end
    end
end
labels_dict = cat(1, labels, labels_flipped)';

% flip if stimulated on the left
fprintf('flipping:\nsubject')
for d = 1:length(dataset.original)
    if contains(dataset.original(d).header.name, 'left')
        fprintf(' %d -', params.subjects(d))
        header = dataset.original(d).header;
        data = dataset.original(d).data;
        [dataset.flipped(d).header, dataset.flipped(d).data, ~] = RLW_flip_electrodes(header, data, labels_dict);
    else
        dataset.flipped(d).header = dataset.original(d).header;
        dataset.flipped(d).data = dataset.original(d).data;
    end
end
clear d i labels labels_flipped labels_dict electrode_n label_new data header

% average across trials and export in ASCII format
fields = fieldnames(dataset); 
for f = 1:length(fields)
    fprintf('exporting %s data:\nsubject ', fields{f})

    % create a separate folder for export, if necessary
    if exist(sprintf('%s\\export\\%s', folder.output, fields{f})) == 0
        mkdir(sprintf('%s\\export', folder.output), fields{f})
    end

    % get the subset
    statement = sprintf('subset = dataset.%s;', fields{f});
    eval(statement)

    % export individual averages
    for d = 1:length(subset)
        % provide update
        fprintf('%d ...', d)
    
        % determine segment name
        if params.subjects(d) < 10
            subject_ID = sprintf('S00%d', params.subjects(d));
        else
            subject_ID = sprintf('S0%d', params.subjects(d));
        end
        if f == 1
            segment_name = [subject_ID '_' params.side{d}];
        elseif f == 2
            if strcmp(params.side{d}, 'right')
                segment_name = [subject_ID '_right'];
            elseif strcmp(params.side{d}, 'left')
                segment_name = [subject_ID '_right_flipped'];
            end
        end
    
        % define export name
        export_name = sprintf('%s\\export\\%s\\%s.avr', folder.output, fields{f}, segment_name);
    
        % define export headers
        export_line{1} = sprintf('Npts= %d TSB= %d DI= %d SB= 1.000 SC= 200.0 Nchan= %d SegmentName= %s', ...
            size(subset(d).data, 6), subset(d).header.xstart*1000,  subset(d).header.xstep*1000, size(subset(d).data, 2), segment_name);
        export_line{2} = strjoin({subset(d).header.chanlocs.labels}, ' ');
    
        % average trials
        data = squeeze(mean(subset(d).data, 1));
    
        % write into an ASCII file
        fileID = fopen(export_name, 'w');
        fprintf(fileID, '%s\n', export_line{1});  
        fprintf(fileID, '%s\n', export_line{2});  
        for e = 1:size(data, 1)
            fprintf(fileID, '%f ', data(e, :));  
            fprintf(fileID, '\n');  
        end
        fclose(fileID); 
    end
    fprintf('\n')
end
fprintf('done.\n\n')
clear d e f fields statement subset subject_ID export_name segment_name export_line fileID data


clear params 
