function results_create_figures(options,layers,settings,SPMmat,varargin)
% create_figures(options,layers,settings,SPMmat [,drugcodes])
% 
% DESCRIPTION
% Runs slice_display toolbox to create result figures with general settings
% specified in layers and settings structs. Specific settings are added to
% these structs in this function.
% 
% Whether MRI or PET images are being created is determined based on
% whether there is a directory called "mri" or "pet" (case-insensitive) in
% the path to the SPM.mat. For MRI images a distinction is made between
% individual level and group level images for the figure titles.
%
% Creating the figures requires binary nifti images that contain the
% significant voxels at the specified significance levels. If these do not
% exist, they are created.
% 
% INPUT
% - options   : struct with user specified info
% - layers    : struct with layer info for slice_display toolbox
% - settings  : struct with settings for slice_display toolbox
% - SPMmat    : char; path to SPM.mat
% Optional:
% - drugcodes : dataset with drug unblinding codes
% 
% OUTPUT
% Saved png images of results figures in </path/to/dir/of/SPM>/figures.
% 
% -------------------------------------------------------------------------
% Ruben van den Bosch
% Donders Institute, Radboud University Nijmegen
% July 2019
% 

% Get drugcodes from input, if provided
if ~isempty(varargin)
    drugcodes = varargin{1};
end

% Directories
% -------------------------------------------------------------------------
[dirs.spm,~,~] = fileparts(SPMmat);
dirs.figures   = fullfile(dirs.spm,'figures');

% Create figures dir if necessary
if ~exist(dirs.figures,'dir')
    mkdir(dirs.figures)
end

% Determine modality
if contains(dirs.spm, [filesep 'mri' filesep], 'IgnoreCase',true)
    modality = 'mri';
elseif contains(dirs.spm, [filesep 'pet' filesep], 'IgnoreCase',true)
    modality = 'pet';
end

% Get p threshold as string to use in names
p = regexp(num2str(options.todo.significance.threshold), '\.', 'split');
p = p{2};

% Get all contrast names from options struct
cons = fieldnames(options.todo.contrast);

% Load SPM for contrast info
load(SPMmat);

% Loop over contrasts and process those selected in options.todo
% -------------------------------------------------------------------------
for iCon = 1:numel(SPM.xCon)

    % Get contrast name
    contrastName = SPM.xCon(iCon).name;
    
    % Replace potential white spaces in contrast name with '_'
    contrast = strrep(contrastName, ' ', '_');
    
    % Skip this contrast if not selected
    % .....................................................................
    % There may be more contrasts in the SPM than the main one, e.g. for
    % covariates or negative contrasts. 
    % Check whether to process current contrast based on the main contrast
    % name. Names of covariate contrasts etc must thus include the main
    % contrast name.
    for i = 1:numel(cons)
        if contains(contrast,cons{i})
            % Get contrast basename
            mainConName = char(regexp(contrast,cons{i},'match'));
            break
        end
    end
    % Skip contrast if main contrast was not selected for processing
    if ~options.todo.contrast.(mainConName).do
        continue
    end
        
    % Contrast specific layer settings
    % ---------------------------------------------------------------------
    if strcmp(options.todo.figType,'Tmap')

        % T-map image
        layers(2).color.file = fullfile(dirs.spm,sprintf('spmT_%.4d.nii',iCon));

        % Binary image of significant voxels.
        % Try combined positive and negative contrast significant voxels
        layers(2).mask.file  = fullfile(dirs.spm,sprintf('significant_voxels_%s_%s_p%s.nii', ...
                                                            contrast,options.todo.significance.thresholdType,p));

        % If the binary image file does not exist, create it
        if ~exist(layers(2).mask.file,'file')
            create_significant_voxels_binary(SPMmat,cellstr(contrast),modality,options);
        end
    elseif strcmp(options.todo.figType,'dualCoded')
           
        % Color file and label (escape '_' characters in label to prevent
        % unwanted subscripts)
        layers(2).color.file  = fullfile(dirs.spm,sprintf('con_%.4d.nii',iCon));
        layers(2).color.label = sprintf('beta_{%s} (a.u.)',strrep(contrast,'_','\_'));

        % Opacity file
        layers(2).opacity.file = fullfile(dirs.spm,sprintf('spmT_%.4d.nii',iCon));

        % Layer 3: Contour of significantly activated voxels
        % Try combined positive and negative contrast significant voxels
        layers(3).color.file = fullfile(dirs.spm,sprintf('significant_voxels_combinedPosNeg_%s_%s_p%s.nii', ...
                                            contrast,options.todo.significance.thresholdType,p));
        
        % If that doesn't exist, use pos/neg only contrast binary
        if ~exist(layers(3).color.file,'file')
            layers(3).color.file = fullfile(dirs.spm,sprintf('significant_voxels_%s_%s_p%s.nii', ...
                                                              contrast,options.todo.significance.thresholdType,p));
        end
        
        % If that binary image file does not exist, create it
        if ~exist(layers(3).color.file,'file')
            create_significant_voxels_binary(SPMmat,cellstr(contrast),modality,options);
        end
    end

    % Create and save figures
    % ---------------------------------------------------------------------
    % Get orientations and slices for this contrast.
    orientations = options.todo.contrast.(mainConName).orientations;
    slices       = options.todo.contrast.(mainConName).slices;
    
    % Create figure title
    % .....................................................................
    % If MRI images, distinguish between individual and group level; if PET
    % images use only contrast name.
    switch modality
        case 'mri'
            % Individual level
            if contains(dirs.spm,options.io.firstLevelDir)

                % Get subject and session
                subject = regexp(dirs.spm,'sub-\d\d\d','match');
                subNr   = str2double(regexp(subject{1},'\d\d\d','match'));
                session = regexp(dirs.spm,'ses-drug\d','match');
                sesNr   = str2double(regexp(session{1},'\d','match'));

                % Get drug for this session
                drug = drugcodes.(sprintf('session%d',sesNr))(drugcodes.subject==subNr);

                % Figure title. Escape '_' to prevent subscripts
                settings.fig_specs.title = [subject{1} ' ' session{1} ' ' sprintf('(%s)',drug{1}) ' : ' strrep(contrastName,'_','\_')];

            % Group level. Escape '_' to prevent subscripts
            elseif contains(dirs.spm,options.io.groupLevelDir)

                % Get drug or drug comparison name, or sess average
                [~,dirNm,~] = fileparts(dirs.spm);
                prefix      = regexp(dirNm,'_','split');
                prefix      = prefix{1};

                % Use original name not including underscores for covariates
                settings.fig_specs.title = ['groupLevel : ' prefix ' ' strrep(contrastName,'_','\_')];
            end
        case 'pet'
            % Escape '_' to prevent subscripts
            settings.fig_specs.title = strrep(contrastName,'_','\_');
    end
    
    % Loop over requested orientations
    for iOr = 1:numel(orientations)

        % Orientation specific settings
        settings.slice.orientation              = orientations{iOr};
        settings.slice.disp_slices              = slices.(orientations{iOr});
        if isempty(options.figure.num_columns)
            settings.fig_specs.n.slice_column   = numel(slices.(orientations{iOr}));
        end

        % Create and save figure
        % .................................................................
        % Prevent warnings spam
        warning off
        [~,~,h_figure] = sd_display(layers,settings);

        % Prevent white text becoming black after save
        set(h_figure,'InvertHardCopy','off')
        
        % Save and close figure
        fname = sprintf('result_%s_%s_%s_%s_p%s.png', ...
                            options.todo.figType,contrast,orientations{iOr},options.todo.significance.thresholdType,p);
        saveas(h_figure,fullfile(dirs.figures,fname));
        close(h_figure);
        warning on
    end
end
end