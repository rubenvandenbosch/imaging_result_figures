function create_significant_voxels_binary(SPMmat,contrast,options)
% function to export binary images of significant voxels
% 
% INPUT
% SPMmat    : char; path to SPM.mat file
% contrast  : cellstr; name of contrasts to export significant voxels for.
%             Can be one contrast or an array to export multiple.
% options   : struct with fields:
%               todo.significance.thresholdType : 'uncorrected' or 'fwe'
%               todo.significance.threshold     : p-value threshold
%               todo.significance.extent        : minumum cluster size
% 
% OUTPUT
% Binary nifti files per contrast in which only the significant voxels at
% the specified (in options) significance level are set to 1.
% 
% -------------------------------------------------------------------------
% Ruben van den Bosch
% Donders Institute, Radboud University Nijmegen
% September 2019
%

% Export binary image of significant voxels in contrasts
% =========================================================================
% Create SPM job
% -------------------------------------------------------------------------
clear jobs

% Load SPM for contrast info
load(SPMmat)

% Fill in info on selected contrasts
% .........................................................................
% Loop over existing contrasts in the SPM
conCount = 1;
for iCon = 1:numel(SPM.xCon)

    % Skip this contrast if not selected for export
    if ~(strcmp(SPM.xCon(iCon).name,contrast{conCount}))
        continue
    end

    % Path to SPM
    jobs{conCount}.spm.stats.results.spmmat = cellstr(SPMmat);

    % Empty title string
    jobs{conCount}.spm.stats.results.conspec(1).titlestr = '';

    % Contrast index
    jobs{conCount}.spm.stats.results.conspec(1).contrasts = iCon;

    % Threshold type, value and min cluster size
    if strcmpi(options.todo.significance.thresholdType,'uncorrected')
        jobs{conCount}.spm.stats.results.conspec(1).threshdesc = 'none';
    elseif strcmpi(options.todo.significance.thresholdType,'fwe')
        jobs{conCount}.spm.stats.results.conspec(1).threshdesc = 'fwe';
    end
    jobs{conCount}.spm.stats.results.conspec(1).thresh = options.todo.significance.threshold;
    jobs{conCount}.spm.stats.results.conspec(1).extent = options.todo.significance.extent;

    % Don't use mask
    jobs{conCount}.spm.stats.results.conspec(1).mask.none = 1;
    
    % Define basename for output image
    % Get p as string for use in basename
    p = regexp(num2str(options.todo.significance.threshold), '\.', 'split');
    p = p{2};
    jobs{conCount}.spm.stats.results.export{1}.binary.basename = sprintf('significant_voxels_%s_%s_p%s', ...
                                                                            SPM.xCon(iCon).name, ...
                                                                            options.todo.significance.thresholdType, ...
                                                                            p);

    % Export as binary
    jobs{conCount}.spm.stats.results.units = 1;
    
    % Increase contrast count
    conCount = conCount + 1;
    
    % Break if all requested contrasts have been processed
    if conCount > numel(contrast)
        break
    end
end

% Run SPM job
% -------------------------------------------------------------------------
jobName = 'export_contrasts_binary';
run_spm_jobs(jobName,jobs)

% SPM prepends "spmT_xxxx" to exported files. Remove it from the binary
% output files, because they are not Tmaps.
% -------------------------------------------------------------------------
[outDir,~,~] = fileparts(SPMmat);
files = dir(fullfile(outDir,'spmT_*significant_voxels*'));
for ifile = 1:numel(files)
    newNm = regexprep(files(ifile).name,'spmT_\d\d\d\d_','');
    movefile(fullfile(files(ifile).folder,files(ifile).name), fullfile(files(ifile).folder,newNm));
end
end

function run_spm_jobs(jobName,jobs)
% RUN_SPM_JOBS(jobName,jobs)
% 
% Subfunction to run spm jobs.
% 

% Initialise cgf utilities
spm_jobman('initcfg');

% Initialise job
jobId = cfg_util('initjob', jobs);

% If successful run job
sts = cfg_util('isjob_id', jobId);
if sts
    cfg_util('run', jobId);
else
    error('Error in initialising %s job.',jobName)
end
end