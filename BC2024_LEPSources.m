%% BC2024: Cortical sources of laser-evoked potentials
% ------------------------------------------------------------------------
% author:   Dominika Sulcova
%           MSH - Medical School Hamburg
% created:  November 2024   
% student:  Pauline Kraus
% ------------------------------------------------------------------------
% project:  
% 
% data:     
% 
% script:   
%           
% output:   
% 

%% params
% directories
folder.LEP = uigetdir(pwd, 'Choose the folder with LEP data');          % MSH D-drive folder --> processed LEPs 
folder.output = uigetdir(pwd, 'Choose the output folder');              % output folder --> local folder with figures, output files, export files
cd(folder.output)

% output
study = 'LEPSources';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

%% select participants 
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
% load single-trial data
dataset = struct;
for s = 1:length(params.subjects)
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
    dataset(s).header = header;
    dataset(s).data = data;   
end

% average across trials and export in ASCII format
fprintf('exporting:\nsubject ')
for d = 1:length(dataset)
    % update
    fprintf('%d ...', params.subjects(d))

    % determine segment name
    if params.subjects(d) < 10
        segment_name = sprintf('S00%d_%s', params.subjects(d), params.side{d});
    else
        segment_name = sprintf('S0%d_%s', params.subjects(d), params.side{d});
    end

    % define export name
    export_name = sprintf('%s\\export\\%s.avr', folder.output, segment_name);

    % define export headers
    export_line{1} = sprintf('Npts= %d TSB= %d DI= %d SB= 1.000 SC= 200.0 Nchan= %d SegmentName= %s', ...
        size(dataset(d).data, 6), dataset(d).header.xstart*1000,  dataset(d).header.xstep*1000, size(dataset(d).data, 2), segment_name);
    export_line{2} = strjoin({dataset(d).header.chanlocs.labels}, ' ');

    % average trials
    data = squeeze(mean(dataset(d).data, 1));

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
fprintf('done.\n')


clear params b d e s subject_ID dataset_block data2load data header export_name segment_name export_line fileID

