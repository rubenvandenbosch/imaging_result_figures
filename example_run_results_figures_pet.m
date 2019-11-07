% example_run_results_figures_pet
% 
% Example script to create and save results figures for a voxelwise PET-BEH
% analysis. Adjust input-output paths and other options in first part of 
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
                  'derivBEHdir', fullfile(projectDir,'bids','derivatives','pet'), ...
                  'derivPETBEHdir', fullfile(projectDir,'bids','derivatives','beh','task','groupLevel','pet'), ...
                  'codeDir', fullfile(projectDir,'code','analysis'), ...
				  'functionsDir', fullfile(projectDir,'code','analysis','mri'), ...
				  'slice_display_dir', fullfile(projectDir,'code','analysis','slice_display'), ...
                  'panel_dir', fullfile(projectDir,'code','analysis','panel'), ...
                  'codeDirMRI', fullfile(projectDir,'code','analysis','mri'));

% Steps to perform
% -------------------------------------------------------------------------
todo.createFigures  = true;
todo.createHTML     = true;

% Output figure type: Tmap, dualCoded
% -------------------------------------------------------------------------
todo.figType = 'dualcoded';

% Covariate settings
% -------------------------------------------------------------------------
% Name covariate "contrast" in options struct, for use in
% results_create_figures function.
% Specify which covariates to make figures for and select orientations and
% which slices to display for each covariate. For options that are not
% specified the defaults (in todo.defaults) are used.
% 
% todo.contrast.<contrastName>.do               : true/false
% todo.contrast.<contrastName>.orientations     : cell array of strings
% todo.contrast.<contrastName>.slices.axial     : vector of slice numbers
% todo.contrast.<contrastName>.slices.coronal   : vector of slice numbers
% todo.contrast.<contrastName>.slices.sagittal  : vector of slice numbers
% todo.contrast.<contrastName>.orientations_in_html : array of cellstr

% Defaults based on an FDOPA study, so focused on striatum
todo.defaults.orientations         = {'axial','coronal','sagittal'};
todo.defaults.slices.axial         = -8:4:28;
todo.defaults.slices.coronal       = -20:4:26;
todo.defaults.slices.sagittal      = [-32:4:-8 8:4:34];
todo.defaults.orientations_in_html = {'axial','coronal','sagittal'};

todo.contrast.<contrast_name>.do        = true;

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

% Optionally set a maximum number of columns for the output figures. Leave
% empty, [], for unspecified.
figure.num_columns = 5;

% Combine all in options struct
% =========================================================================
options.io        = iostruct;
options.todo      = todo;
options.figure    = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get figure type name with correct capitals for use in filenames
if strcmpi(options.todo.figType,'tmap')
    options.todo.figType = 'Tmap';
elseif strcmpi(options.todo.figType,'dualcoded')
    options.todo.figType = 'dualCoded';
end

% Add path with functions for figure and html creation
addpath(options.io.functionsDir)

% =========================================================================
% CREATE FIGURES
% =========================================================================

if options.todo.createFigures
    
    % Make sure directory containing slice display functions is on path
    sd_dir = fullfile(options.io.codeDir,'slice_display');
    addpath(genpath(sd_dir));

    % Get custom colormaps
    load(fullfile(sd_dir,'colormaps.mat'));

    % Define common settings in layers
    % ---------------------------------------------------------------------
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

        if strcmpi(options.todo.significance.thresholdType,'uncorrected') && options.todo.significance.threshold == 0.001
            layers(2).opacity.range     = [0 3.14];
        elseif strcmpi(options.todo.significance.thresholdType,'uncorrected') && options.todo.significance.threshold == 0.01
            layers(2).opacity.range     = [0 2.33];
        elseif strcmpi(options.todo.significance.thresholdType,'fwe') && options.todo.significance.threshold == 0.05
            layers(2).opacity.range     = [0 3.18];
        end
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
    
    % Loop over contrasts
    % ---------------------------------------------------------------------
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
        
        % Path to SPM
        SPMmat = fullfile(options.io.derivPETBEHdir,conName,'SPM.mat');
        
        % Create figures
        results_create_figures(options,layers,settings,SPMmat)
        
    end
end

% =========================================================================
% CREATE HTML
% =========================================================================

if options.todo.createHTML
    
    % Get p threshold as string to use in names
    p = regexp(num2str(options.todo.significance.threshold), '\.', 'split');
    p = p{2};
    
    % Create HTML per contrast with the three (or available) orientations
    % and create HTML with collection of all contrast results
    % ---------------------------------------------------------------------
    collectionFiles = {};
    
    % Loop over contrasts
    cons = fieldnames(options.todo.contrast);
    for icon = 1:numel(cons)
        conName = cons{icon};
    
        % Skip contrast if not selected in options.todo
        if ~options.todo.contrast.(conName).do
            continue
        end

        % Create HTML with three orientations for each contrast separately
        % .................................................................
        % HTML output file
        outputDir = fullfile(options.io.derivPETBEHdir,conName,'figures');
        if ~exist(outputDir,'dir')
            continue
        end
        
        htmlFile = fullfile(outputDir,sprintf('result_%s_%s_%s_p%s.html',...
                                                options.todo.figType,conName,options.todo.significance.thresholdType,p));

        % Collect only the regular, not inverse, images to put in html
        imgFiles = cellstr(spm_select('FPList',outputDir,sprintf('^result_%s_%s_.*_%s_p%s.png$', ...
                                                                    options.todo.figType,conName,options.todo.significance.thresholdType,p)));
        
        % Title and headings
        title.title = sprintf('%s: PET-BEH results',conName);
        title.head1 = sprintf('Voxelwise PET-BEH results for %s',conName);
        title.head2 = sprintf('Significant voxels at: p = %s, %s',num2str(options.todo.significance.threshold),options.todo.significance.thresholdType);

        % Search string for image file paths to base separation of images 
        % with big and small borders on. 
        % Set to false for no border
        borders.big     = false;
        borders.small   = false;

        % Create HTML
        results_create_html(htmlFile,title,imgFiles,borders)
        
        % Add image files to the list for the all contrasts collection HTML
        % .................................................................
        collectionFiles = [collectionFiles; imgFiles];
    end
    
    % Create HTML with collection of all contrast results
    % ---------------------------------------------------------------------
    % HTML file
    outputDir = fullfile(options.io.derivPETBEHdir,'figure_collections');
    if ~exist(outputDir,'dir')
        mkdir(outputDir)
    end
    htmlFile = fullfile(outputDir,sprintf('results_%s_allCovariates_%s_p%s.html',...
                                            options.todo.figType,options.todo.significance.thresholdType,p));
    
    % Title and headings
    title.title = 'Voxelwise PET-BEH results';
    title.head1 = 'Voxelwise PET-BEH results: collection of all covariates';
    title.head2 = sprintf('Significant voxels at: p = %s, %s',num2str(options.todo.significance.threshold),options.todo.significance.thresholdType);

    % Search string for image file paths to base separation of images 
    % with big and small borders on. 
    % Set to false for no border
    borders.big     = 'pCorrect';
    borders.small   = false;

    % Create HTML
    results_create_html(htmlFile,title,collectionFiles,borders)
end
