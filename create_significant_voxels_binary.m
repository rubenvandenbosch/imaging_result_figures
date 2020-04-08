function create_significant_voxels_binary(SPMmat,contrast,modality,sigOptions,varargin)
% create_significant_voxels_binary(SPMmat,contrast,modality,sigOptions)
% 
% SUMMARY
% Function to export binary images of significant voxels. The exported
% binary images include the contrast's significant voxels in both postive 
% and negative direction. Because SPM only calculates statistics on effects
% of a contrast in the positive direction, the inverse of the contrast is
% created and it's significant voxels are combined with the positive
% direction's significant voxels to create the image with significant
% voxels in both directions.
% 
% INPUT
% SPMmat     : char; path to SPM.mat file
% contrast   : cellstr; name of contrasts to export significant voxels for.
%              Can be one contrast or an array to export multiple.
% modality   : char; 'mri' OR 'pet'
% sigOptions : struct with at least these required fields:
%       thresholdType : char; 'uncorrected' OR 'fwe'
%       threshold     : cell array with p-value threshold(s)
%       extent        : double; minumum cluster size
% 
% Optional input:
% create_significant_voxels_binary(SPMmat,contrast,modality,sigOptions,isInitCfg)
% 
% isInitCfg  : true/false; have SPM cfg utilities already been loaded. If
%              true they are not loaded again to save time. If not
%              specified the utilities will be loaded.
% 
% OUTPUT
% Binary nifti files per contrast in which only the significant voxels (for
% both positive and negative effects) at the requested significance level 
% (specified in sigOptions) are set to 1.
% 
% -------------------------------------------------------------------------
% Ruben van den Bosch
% Donders Institute, Radboud University Nijmegen
% September 2019
%

% Process optional input argument
if ~isempty(varargin)
    assert(numel(varargin) == 1, 'Too many input arguments');
    isInitCfg = varargin{1};
else
    isInitCfg = false;
end

% Create SPM jobs
% =========================================================================
clear jobs

% Load modality defaults
if strcmpi(modality,'MRI')
    modality = 'FMRI';
end
spm('defaults',upper(modality));

% Initialise SPM cfg utilities if necessary
if ~isInitCfg
    spm_jobman('initcfg');
end

% Load SPM and get contrasts info
% -------------------------------------------------------------------------
load(SPMmat)

% Get contrast names of all existing contrasts in the SPM
conNms = {SPM.xCon(:).name}';

% Replace potential white spaces in SPM contrast names with "_"
conNms = strrep(conNms,' ','_');

% Replace potential white spaces in contrast names requested for processing
% with "_". Prevents white space differences in SPM contrast names and
% requestion names from messing up name comparisons etc.
contrast = strrep(contrast,' ','_');

% Loop over existing contrasts in the SPM and process only selected
% contrasts
% -------------------------------------------------------------------------
conCount = 1;  % counter for requested contrasts processed
jobCount = 1;  % counter for job numbers

for iCon = 1:numel(conNms)
    conName = conNms{iCon};
    
    % Skip this contrast if not selected for export or if it's a negative
    % contrast of another requested contrast
    if ~ismember(cellstr(conName),contrast)
        continue
    elseif startsWith(conName,'negative_') && ismember(cellstr(strrep(conName,'negative_','')),contrast)
        continue
    elseif startsWith(conName,'inverse_') && ismember(cellstr(strrep(conName,'inverse_','')),contrast)
        continue
    end
    
    % Check whether the negative, i.e. inverse, contrast exists. If not
    % create that negative contrast, so we can later combine it's
    % significant voxels with those of the regular contrast to have both
    % directions in the binary image.
    % ---------------------------------------------------------------------
    % The negative, inverse, contrast would have "negative_" or "inverse_"
    % as prefix to the contrast name
    if ~ismember(['negative_' conName],conNms) && ~ismember(['inverse_' conName],conNms)
        
        % Make job to create negative contrast
        % .................................................................
        clear newConJob
        newConJob{1}.spm.stats.con.spmmat = cellstr(SPMmat);
        
        % Contrast name. Use original contrast names, with white spaces
        % Set negative prefix for contrast
        negPrefix = 'negative_';
        newConJob{1}.spm.stats.con.consess{1}.tcon.name = [negPrefix SPM.xCon(iCon).name];
        
        % New contrast weights
        newConJob{1}.spm.stats.con.consess{1}.tcon.weights  = -SPM.xCon(iCon).c;
        
        % No repeat over sessions
        newConJob{1}.spm.stats.con.consess{1}.tcon.sessrep  = 'none';
        
        % Do not delete existing contrasts
        newConJob{1}.spm.stats.con.delete = 0;
        
        % Run job to create negative contrast. Must be run before other
        % jobs, because this new contrast's index number is determined from
        % the SPM struct later for which it has to have been added to SPM.
        % .................................................................
        jobName = 'create_negative_contrast';
        run_spm_jobs(jobName,newConJob)
        
        % After running the contrast creation job, reload SPM.mat to have
        % the updated SPM struct in working memory
        load(SPMmat);
        
    % If the negative contrast exists, but with the prefix "inverse_" then 
    % no need to create new contrast, but store prefix for use below
    elseif ismember(['negative_' conName],conNms)
        negPrefix = 'negative_';
    elseif ismember(['inverse_' conName],conNms)
        negPrefix = 'inverse_';
    end
    
    % Create jobs to write both the positive and negative contrasts'
    % significant voxels to binary files
    % ---------------------------------------------------------------------
    % Loop over potential multiple requested p-thresholds for export
    if ~iscell(sigOptions.threshold)
        sigOptions.threshold = {sigOptions.threshold};
    end
    for ip = 1:numel(sigOptions.threshold)

        % Get p-value threshold as string for use in file names
        p = regexp(num2str(sigOptions.threshold{ip}), '\.', 'split');
        p = p{2};
        
        % Loop twice to create binary images for both positive and negative
        % effects
        for ix = 1:2

            % Get contrast name with correct prefix (for negative one), and
            % set the outname to positive or negative prefix.
            if ix == 1
                name    = conName;
                outname = ['positive_' conName];
            elseif ix == 2
                name    = [negPrefix conName];

                % For file output always use "negative_" prefix instead of
                % "inverse_"
                outname = ['negative_' conName];
            end

            % Path to SPM
            jobs{jobCount}.spm.stats.results.spmmat = cellstr(SPMmat);

            % Empty title string
            jobs{jobCount}.spm.stats.results.conspec(1).titlestr = '';

            % Find contrast index. Use names from SPM struct because new
            % contrasts may have been added in previous jobs. Do replace
            % potential white spaces in SPM struct contrast names with "_"
            % when comparing with current contrast name.
            for conIx = 1:numel(SPM.xCon)
                if strcmp(name,strrep(SPM.xCon(conIx).name,' ','_'))
                    jobs{jobCount}.spm.stats.results.conspec(1).contrasts = conIx;
                    break
                end
            end

            % Threshold type, value and min cluster size
            if strcmpi(sigOptions.thresholdType,'uncorrected')
                jobs{jobCount}.spm.stats.results.conspec(1).threshdesc = 'none';
            elseif strcmpi(sigOptions.thresholdType,'fwe')
                jobs{jobCount}.spm.stats.results.conspec(1).threshdesc = 'fwe';
            end
            jobs{jobCount}.spm.stats.results.conspec(1).thresh = sigOptions.threshold{ip};
            jobs{jobCount}.spm.stats.results.conspec(1).extent = sigOptions.extent;

            % Don't use mask
            jobs{jobCount}.spm.stats.results.conspec(1).mask.none = 1;

            % Define basename for output image.
            basename = sprintf('significant_voxels_%s_%s_p%s',outname,sigOptions.thresholdType,p);
            jobs{jobCount}.spm.stats.results.export{1}.binary.basename = basename;

            % Store the fullpath filename of written image to use below in 
            % the combination job. SPM automatically adds prefix "spmT_" to
            % the written files.
            [outDir,~,~] = fileparts(SPMmat);
            if ix == 1
                pos = fullfile(outDir,sprintf('spmT_%.4d_%s.nii',conIx,basename));
            elseif ix ==2
                neg = fullfile(outDir,sprintf('spmT_%.4d_%s.nii',conIx,basename));
            end

            % Export as binary
            jobs{jobCount}.spm.stats.results.units = 1;

            % Increase jobCount
            jobCount = jobCount + 1;
        end

        % Create job to combine postive and negative direction significant
        % voxels into one binary image
        % -----------------------------------------------------------------
        % Output dir
        [outDir,~,~] = fileparts(SPMmat);

        % Output file name
        outputname = sprintf('significant_voxels_%s_%s_p%s',conName,sigOptions.thresholdType,p);

        % Fill imcalc job
        % Input files are the pos and neg files; names stored in loop above
        jobs{jobCount}.spm.util.imcalc.input            = cellstr(char(pos,neg));
        jobs{jobCount}.spm.util.imcalc.output           = outputname;
        jobs{jobCount}.spm.util.imcalc.outdir           = cellstr(outDir);
        jobs{jobCount}.spm.util.imcalc.expression       = 'i1 | i2';
        jobs{jobCount}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
        jobs{jobCount}.spm.util.imcalc.options.dmtx     = 0;
        jobs{jobCount}.spm.util.imcalc.options.mask     = 0;
        jobs{jobCount}.spm.util.imcalc.options.interp   = 0;
        jobs{jobCount}.spm.util.imcalc.options.dtype    = 2;

        % Increase jobCount
        jobCount = jobCount + 1;
    end
    
    % Increase counter of requested contrasts processed
    conCount = conCount + 1;
    
    % Break if all requested contrasts have been processed
    if conCount > numel(contrast)
        break
    end
end

% Run SPM jobs
% =========================================================================
jobName = 'export_contrasts_binary';
run_spm_jobs(jobName,jobs)

% Remove the separate positive and negative significant voxels images
% =========================================================================
[outDir,~,~] = fileparts(SPMmat);
files = cellstr(char(spm_select('FPList',outDir,'^spmT_\d\d\d\d_significant_voxels_positive_.*.nii$'), ...
                     spm_select('FPList',outDir,'^spmT_\d\d\d\d_significant_voxels_negative_.*.nii$')));
for ifile = 1:numel(files)
    delete(files{ifile});
end

% SPM prepends "spmT_xxxx" to exported files. Remove it from the binary
% output files, because they are not Tmaps.
% =========================================================================
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