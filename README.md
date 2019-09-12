# imaging_result_figures

## Dependencies
- [MATLAB R2015a or later](http://www.mathworks.com)
- [MATLAB Statistics Toolbox](http://www.mathworks.com/products/statistics/)
- [SPM](http://www.fil.ion.ucl.ac.uk/spm/) by the Wellcome Trust Centre for Neuroimaging at UCL
- [slice_display toolbox](https://github.com/bramzandbelt/slice_display) by Bram Zandbelt
- [Panel](https://www.mathworks.com/matlabcentral/fileexchange/20003-panel) by Ben Mitch

**Note:** Adapt the function sd_display (from slice_display toolbox) line 1 to: `function [settings, p, h_figure] = sd_display(layers, settings)`. A pull request has been made to slice_display for this change.

## Description
Code to create and save result figures from task fMRI analyses or PET-BEH voxelwise correlations.

Requires first level or group level analyses to have been run using SPM.

Adapt the user input in the example code to your needs in order to run it.

## Outputs
- png images per orientation in each contrast's directory
- html collection of the images in different orientations in each contrast's directory
- html collection of all figures in a separate directory (can specify which orientations to include)

## Author
Ruben van den Bosch
Donders Institute, Radboud University Nijmegen
The Netherlands
