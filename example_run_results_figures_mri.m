% example_run_results_figures_mri
% 
% Example script to create and save results figures for a task-fMRI 
% experiment. Adjust input-output paths and other options in first part of 
% the file, until "END USER INPUT"
% 

clear

% Directories
% -------------------------------------------------------------------------
if isunix
    projectDir = '</path/to/project/directory>';
    addpath </path/to/spm12>
elseif ispc
    projectDir = '<\path\to\projectdirectory>';
    addpath <\path\to\spm12>
end

% Example of input output directory structure
iostruct = struct('projectDir', projectDir, ...
                  'bidsDir', fullfile(projectDir,'bids'), ...
                  'derivDir', fullfile(projectDir,'bids','derivatives'), ...
                  'derivMRIdir',fullfile(projectDir,'bids','derivatives','mri'), ...
                  'firstLevelDir',fullfile(projectDir,'bids','derivatives','mri','firstLevel'), ...
                  'groupLevelDir',fullfile(projectDir,'bids','derivatives','mri','groupLevel'), ...
                  'codeDir', fullfile(projectDir,'code','analysis'), ...
				  'functionsDir', fullfile(projectDir,'code','analysis','mri'), ...
                  'slice_display_dir', fullfile(projectDir,'code','analysis','slice_display'), ...
                  'panel_dir', fullfile(projectDir,'code','analysis','panel'), ...
                  'codeDirFirstLevel', fullfile(projectDir,'code','analysis','mri','1stLevel'), ...
                  'codeDirGroupLevel', fullfile(projectDir,'code','analysis','mri','2ndLevel'));

% Individual level or group level?
% -------------------------------------------------------------------------
% Cellstr of levels to process; {'individual','group'}
level = {'individual','group'};

% Steps to perform
% -------------------------------------------------------------------------
todo.createFigures  = true;
todo.createHTML     = true;

% Output figure type: Tmap, dualCoded
% -------------------------------------------------------------------------
todo.figType = 'dualCoded';

% Contrast settings
% -------------------------------------------------------------------------
% Specify which contrasts to make figures for and select orientations and
% which slices to display for each contrast.
% 
% todo.contrast.<contrastName>.do               : true/false
% todo.contrast.<contrastName>.orientations     : cell array of strings
% todo.contrast.<contrastName>.slices.axial     : vector of slice numbers
% todo.contrast.<contrastName>.slices.coronal   : vector of slice numbers
% todo.contrast.<contrastName>.slices.sagittal  : vector of slice numbers
% todo.contrast.<contrastName>.orientations_in_html : cell array of strings

% Default settings to use if orientation and slice settings have not been
% specified for a contrast that is set to true (here focused on striatum)
todo.defaults.orientations         = {'axial','coronal','sagittal'};
todo.defaults.slices.axial         = -8:4:28;
todo.defaults.slices.coronal       = -20:4:26;
todo.defaults.slices.sagittal      = [-32:4:-8 8:4:34];
todo.defaults.orientations_in_html = {'axial','coronal','sagittal'};

% Example for right-hand buttonPress contrast
todo.contrast.buttonPress.do                    = true;
todo.contrast.buttonPress.orientations          = {'axial','coronal','sagittal'};
todo.contrast.buttonPress.slices.axial          = 30:2:46;
todo.contrast.buttonPress.slices.coronal        = -68:6:-20;
todo.contrast.buttonPress.slices.sagittal       = [-60:4:-34 0 50];
todo.contrast.buttonPress.orientations_in_html  = {'axial','coronal','sagittal'};

% This code was developed for a multi-session drug study. The drug names 
% are included in the output for MRI grouplevel results and can thus be 
% used to select based on those names. Specify drug names and comparisons
% to include here.
% If your setup was different you can make this empty. You'll have to 
% implement the groupLevel file selection yourself in the code below in
% that case.
% -------------------------------------------------------------------------
todo.drugs      = {'drug1','drug2','drug3','drug1vsdrug2','drug1vsdrug3','drug2vsdrug3'};

% Path to drug unblinding codes (first level was coded using session 
% numbers. This is the path to the unblinding table of which session had 
% which drug)
% -------------------------------------------------------------------------
todo.drugcodes = fullfile(projectDir,'</path/to/unblinding/code.csv>');

% Process group level average of sessions? (if multi-session study)
% -------------------------------------------------------------------------
todo.average    = false;

% Significance level
% -------------------------------------------------------------------------
% thresholdType : 'uncorrected', 'fwe'
% threshold     : p-value threshold
% extent        : minimum cluster size in number of voxels to include in
%                 the binary images that are created (if necessary) of the 
%                 significant voxels in a contrast
todo.significance.thresholdType   = 'uncorrected';
todo.significance.threshold       = 0.001;
todo.significance.extent          = 0;

% Which denoise level to use. (ICA-AROMA denoised data, using not 
% nonaggressive, nor aggressive, but the regressors in glm. Possible to 
% include more here to loop over different denoise levels)
% -------------------------------------------------------------------------
todo.denoise.semiaggressive = true;

% Subjects and sessions to include
% -------------------------------------------------------------------------
subjects = 1:100;
sessions = [1 2 3];

% Path to anatomical scan to use to project blobs on
% -------------------------------------------------------------------------
% fig.anatImage = fullfile(spm('Dir'),'canonical','single_subj_T1.nii');
fig.anatImage = fullfile(projectDir,'bids','derivatives','mri','fmriprep','group','group_average_T1scan.nii');

% Set the absolute t-value range that codes opacity in dual-coded images.
% -------------------------------------------------------------------------
% Open a contrast's results window in SPM GUI to see t-threshold.
% Option to set different values for different significance settings.
if strcmpi(todo.significance.thresholdType,'uncorrected') && todo.significance.threshold == 0.001
    fig.opacityRange = [0 3.18];
elseif strcmpi(todo.significance.thresholdType,'uncorrected') && todo.significance.threshold == 0.01
    fig.opacityRange = [0 2.33];
elseif strcmpi(todo.significance.thresholdType,'fwe') && todo.significance.threshold == 0.05
    fig.opacityRange = [0 3.18]; % [0 4.83];
end

% Optionally set a maximum number of columns for the output figures. Leave
% empty, [], for unspecified.
% -------------------------------------------------------------------------
fig.num_columns = [];

% Combine all in options struct
% =========================================================================
options.io        = iostruct;
options.level     = level;
options.todo      = todo;
options.subjects  = subjects;
options.sessions  = sessions;
options.figure    = fig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get figure type name with correct capitals for use in filenames
% -------------------------------------------------------------------------
if strcmpi(options.todo.figType,'tmap')
    options.todo.figType = 'Tmap';
elseif strcmpi(options.todo.figType,'dualcoded')
    options.todo.figType = 'dualCoded';
end

% Add path with functions for figure and html creation
addpath(options.io.functionsDir)

% Fill in missing orientation and slice settings for to be processed
% contrasts with default settings
% -------------------------------------------------------------------------
cons = fieldnames(options.todo.contrast);
for icon = 1:numel(cons)
    conName = cons{icon};

    % Skip this contrast if not selected in options.todo
    if ~options.todo.contrast.(conName).do
        continue
    end

    % Fill in default orientation and slice settings for those settings
    % that were not specified.
    opts = fieldnames(options.todo.defaults);
    for iOpt = 1:numel(opts)
        if ~isfield(options.todo.contrast.(conName),opts{iOpt})
            options.todo.contrast.(conName).(opts{iOpt}) = options.todo.defaults.(opts{iOpt});
        end
    end
end

% =========================================================================
% CREATE FIGURES
% =========================================================================

if options.todo.createFigures
    
    % Loop over denoising levels
    % ---------------------------------------------------------------------
    denoise_lvls = fieldnames(options.todo.denoise);

    for idenoiseLvl = 1:numel(denoise_lvls)

        % Skip if this denoise level is set to false
        if ~options.todo.denoise.(denoise_lvls{idenoiseLvl})
            continue
        end

        % Takes name for level of denoising from options structure
        denoise = denoise_lvls{idenoiseLvl};

        % Add directories with panel and slice display functions to path
        sd_dir = fullfile(options.io.codeDir,'slice_display');
        addpath(genpath(options.io.slice_display_dir));
        addpath(options.io.panel_dir);
        
        % Get custom colormaps
        load(fullfile(sd_dir,'colormaps.mat'));

        % Define common settings in layers
        % -----------------------------------------------------------------
        if strcmp(options.todo.figType,'Tmap')
            layers                          = sd_config_layers('init',{'truecolor','blob'});
        elseif strcmp(options.todo.figType,'dualCoded')
            layers                          = sd_config_layers('init',{'truecolor','dual','contour'});
        end
        settings                            = sd_config_settings('init');

        % Layer 1: anatomical image (same for all, MNI space)
        layers(1).color.file                = fullfile(spm('Dir'),'canonical','single_subj_T1.nii');
        layers(1).color.map                 = gray(256);

        % Layer 2: thresholded t-map or dual-coded
        layers(2).color.map                 = CyBuBkRdYl;

        % Thresholded t-map general options
        if strcmp(options.todo.figType,'Tmap')
            layers(2).color.label           = 't-value';

        % Dual-coded general options
        elseif strcmp(options.todo.figType,'dualCoded')
            layers(2).opacity.label         = '| t |';
            layers(2).opacity.range         = options.figure.opacityRange;
        end

        % Additional range setting for color map to fix/match them
        % layers(2).color.range               = [-0.0069 0.0069];

        % Maximumn number of columns and figure size
        if ~isempty(options.figure.num_columns)
            settings.fig_specs.n.slice_column = options.figure.num_columns;
            if options.figure.num_columns > 7
                settings.fig_specs.width.figure     = 400;
            else
                settings.fig_specs.width.figure     = 300;
            end
        else
            settings.fig_specs.width.figure     = 400;
        end
        
        % Loop over individual or group level todo
        % -----------------------------------------------------------------
        for iLevel = 1:numel(options.level)
            if strcmpi(options.level{iLevel},'individual')
                
                % Load drug codes so drug can be added to figure titles
                % ---------------------------------------------------------
                drugcodes = readtable(options.todo.drugcodes,'Delimiter',',');
                
                % Loop over subjects and sessions
                % ---------------------------------------------------------
                for isub = 1:numel(options.subjects)
                    iSubject = options.subjects(isub);
                    
                    for ises = 1:numel(options.sessions)
                        iSession = options.sessions(ises);

                        % Path to SPM.mat
                        SPMmat = fullfile(options.io.firstLevelDir,sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession),denoise,'SPM.mat');
                        
                        % Skip session if SPM.mat does not exist
                        if ~exist(SPMmat,'file')
                            warning('%s does not exist. Skipping this session.',SPMmat);
                            continue
                        end
                        
                        % Create figures
                        results_create_figures(options,layers,settings,SPMmat,drugcodes)
                    end
                end
                
            elseif strcmpi(options.level{iLevel},'group')
                
                % Contrast names
                cons = fieldnames(options.todo.contrast);
                
                % Loop over contrasts
                % ---------------------------------------------------------
                for icon = 1:numel(cons)
                    
                    contrast = cons{icon};
                    
                    % Skip this contrast if not selected
                    if ~options.todo.contrast.(contrast).do
                        continue
                    end
                    
                    % Loop over drugs and drug comparisons 
                    % .....................................................
                    for idrug = 1:numel(todo.drugs)
                        % Path to group level contrast SPM
                        SPMmat = fullfile(options.io.groupLevelDir,denoise,sprintf('%s_%s',todo.drugs{idrug},contrast),'SPM.mat');

                        % Create figures
                        results_create_figures(options,layers,settings,SPMmat)
                    end
                    
                    % Process sessions average if true
                    % .....................................................
                    if options.todo.average
                        % Path to group level contrast SPM
                        SPMmat = fullfile(options.io.groupLevelDir,denoise,sprintf('sessionAverage_%s',contrast),'SPM.mat');

                        % Skip if SPM.mat does not exist
                        if ~exist(SPMmat,'file')
                            warning('%s does not exist. Skipping this contrast.',SPMmat);
                            continue
                        end
                        
                        % Create figures
                        results_create_figures(options,layers,settings,SPMmat)
                    end
                end
            else
                error("Unrecognized level for analysis. Valid inputs: 'individual' or 'group' or both of those")
            end
        end
    end
end

% =========================================================================
% CREATE HTML
% =========================================================================

if options.todo.createHTML
    
    % Loop over denoising levels
    % ---------------------------------------------------------------------
    denoise_lvls = fieldnames(options.todo.denoise);

    for idenoiseLvl = 1:numel(denoise_lvls)

        % Skip if this denoise level is set to false
        if ~options.todo.denoise.(denoise_lvls{idenoiseLvl})
            continue
        end

        % Takes name for level of denoising from options structure
        denoise = denoise_lvls{idenoiseLvl};
        
        % Get p threshold as string to use in names
        p = regexp(num2str(options.todo.significance.threshold), '\.', 'split');
        p = p{2};

        % Loop over individual or group level todo
        % -----------------------------------------------------------------
        for iLevel = 1:numel(options.level)
            if strcmpi(options.level{iLevel},'individual')

                % Loop over contrasts
                % ---------------------------------------------------------
                cons = fieldnames(options.todo.contrast);
                for iCon = 1:numel(cons)
                    contrast = cons{iCon};

                    % Skip contrast if not selected in options.todo
                    % -----------------------------------------------------
                    if ~options.todo.contrast.(contrast).do
                        continue
                    end
    
                    % HTML output file
                    % -----------------------------------------------------
                    outputDir = fullfile(options.io.firstLevelDir,'figure_collections');

                    % Create outputDir if not exists
                    if ~exist(outputDir,'dir')
                        mkdir(outputDir)
                    end
                    htmlFile = fullfile(outputDir,sprintf('individual_results_%s_%s_%s_%s_p%s.html',...
                                                            options.todo.figType,contrast,denoise,options.todo.significance.thresholdType,p));

                    % Create file list
                    % -----------------------------------------------------
                    imgFiles = {};

                    % Loop over subjects and sessions
                    % .....................................................
                    for isub = 1:numel(options.subjects)
                        iSubject = options.subjects(isub);

                        for ises = 1:numel(options.sessions)
                            iSession = options.sessions(ises);

                            % Subject session directory
                            dirs.firstLevel = fullfile(iostruct.firstLevelDir,sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession),denoise);
                            dirs.figures    = fullfile(dirs.firstLevel,'figures');

                            % Collect image files
                            orientations = options.todo.contrast.(contrast).orientations_in_html;
                            for iOr = 1:numel(orientations)
                                fname = sprintf('result_%s_%s_%s_%s_p%s.png', ...
                                                    options.todo.figType,contrast,orientations{iOr},options.todo.significance.thresholdType,p);
                                imgFiles = [imgFiles; ...
                                            cellstr(fullfile(dirs.figures,fname))];
                            end
                        end
                    end
                    
                    % Title and headings
                    % -----------------------------------------------------
                    title.title = sprintf('%s: individual results',contrast);
                    title.head1 = sprintf('Individual results of contrast: %s',contrast);
                    title.head2 = sprintf('Significant voxels at: p = %s, %s',num2str(options.todo.significance.threshold),options.todo.significance.thresholdType);
                    
                    % Search string for image file paths to base separation
                    % of images with big and small borders on. 
                    % Set to false for no border
                    % -----------------------------------------------------
                    borders.big     = 'sub-';
                    borders.small   = 'ses-drug';
                    
                    % Create HTML with collection of figures
                    % -----------------------------------------------------
                    results_create_html(htmlFile,title,imgFiles,borders)
                end
                
            elseif strcmpi(options.level{iLevel},'group')
                
                % HTML output file
                % ---------------------------------------------------------
                outputDir = fullfile(options.io.groupLevelDir,'figure_collections');

                % Create outputDir is not exists
                if ~exist(outputDir,'dir')
                    mkdir(outputDir)
                end
                
                % Loop over contrasts
                cons = fieldnames(options.todo.contrast);
                for iCon = 1:numel(cons)
                    contrast = cons{iCon};
                    
                    % Skip contrast if not selected in options.todo
                    % -----------------------------------------------------
                    if ~options.todo.contrast.(contrast).do
                        continue
                    end
                    
                    % HTML output files
                    htmlFile = fullfile(outputDir,sprintf('groupLevel_results_%s_%s_%s_%s_p%s.html',...
                                                            options.todo.figType,contrast,denoise,options.todo.significance.thresholdType,p));
                    % Create file list
                    % -----------------------------------------------------
                    imgFiles = {};
                    
                    % Collect image files
                    % .....................................................
                    orientations = options.todo.contrast.(contrast).orientations_in_html;

                    % If orientations to include is empty, skip contrast
                    if isempty(orientations)
                        continue
                    end
                    
                    % Loop over drugs and drug comparisons to include
                    for idrug = 1:numel(options.todo.drugs)
                        for iOr = 1:numel(orientations)
                            fname = sprintf('result_%s_%s_%s_%s_p%s.png', ...
                                                options.todo.figType,contrast,orientations{iOr},options.todo.significance.thresholdType,p);
                            imgFiles = [imgFiles; ...
                                        cellstr(fullfile(options.io.groupLevelDir,denoise,sprintf('%s_%s',options.todo.drugs{idrug},contrast),'figures',fname))];
                        end
                        
                        % Add potential covariate contrasts.
                        % .................................................
                        % Covariate contrasts are named <contrast>_x_<cov>
                        pat = sprintf('result_%s_%s*_%s_%s_p%s.png', ...
                                                options.todo.figType,contrast,orientations{1},options.todo.significance.thresholdType,p);
                        fmatch = dir(fullfile(options.io.groupLevelDir,denoise,sprintf('%s_%s',options.todo.drugs{idrug},contrast),'figures',pat));
                        
                        if numel(fmatch) > 1
                            for icov = 2:numel(fmatch)
                                % Get contrast+covariate name
                                split = strsplit(fmatch(icov).name,'_');
                                conCellNum = find(strcmp(contrast,split));
                                concovNm = [split{conCellNum} '_' split{conCellNum+1} '_' split{conCellNum+2}];
                                
                                % Add files to list
                                for iOr = 1:numel(orientations)
                                    fname = sprintf('result_%s_%s_%s_%s_p%s.png', ...
                                                        options.todo.figType,concovNm,orientations{iOr},options.todo.significance.thresholdType,p);
                                    imgFiles = [imgFiles; ...
                                                cellstr(fullfile(options.io.groupLevelDir,denoise,sprintf('%s_%s',options.todo.drugs{idrug},contrast),'figures',fname))];
                                end
                            end
                        end
                    end
                    
                    % Add sessions average if true
                    % .....................................................
                    if options.todo.average
                        for iOr = 1:numel(orientations)
                            fname = sprintf('result_%s_%s_%s_%s_p%s.png', ...
                                                options.todo.figType,contrast,orientations{iOr},options.todo.significance.thresholdType,p);
                            imgFiles = [imgFiles; ...
                                        cellstr(fullfile(options.io.groupLevelDir,denoise,sprintf('sessionAverage_%s',contrast),'figures',fname))];
                        end
                    end
                
                    % Title and headings
                    % -----------------------------------------------------
                    title.title = sprintf('%s : groupLevel results',contrast);
                    title.head1 = sprintf('GroupLevel results of contrast: %s',contrast);
                    title.head2 = sprintf('Significant voxels at: p = %s, %s',num2str(options.todo.significance.threshold),options.todo.significance.thresholdType);

                    % Search string for image file paths to base separation
                    % of images with big and small borders on. 
                    % Set to false for no border
                    % -----------------------------------------------------
                    borders.big     = contrast;
                    borders.small   = false;
                    
                    % Create HTML with collection of figures
                    % -----------------------------------------------------
                    results_create_html(htmlFile,title,imgFiles,borders)
                end
            end
        end
    end
end