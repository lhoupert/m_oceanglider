# README #

The repository contains the UEA Seaglider data processing toolbox. It is a work in progress, contributions are welcome.

Toolbox branch is the current development branch for the new (as of 2015) version of the toolbox. The old toolbox can still be found under the legacy branch - note that there are numerous bugs and scientific inaccuracies which have been identified but not corrected on the legacy branch.

### How do I get set up? ###

* Summary of set up

1. Pull the "toolbox" branch and add the .m files to your MATLAB search path.

2. data_files = gt_sg_load_merge_data(gliderNumber,'/data/path');

3. data = gt_sg_process_data(data_files)


* Configuration

Glider specific settings are set in sg_calib_constants.m.
General settings are set in gt_sg_settings.m.
gt_sg_settings.m will eventually be provided as a standardised blank template.
Information in sg_calib_constants.m will ALWAYS override the data in gt_sg_settings.m.

* Dependencies

This toolbox requires the Gibbs Seawater toolbox (v.3.04). It must also be added to your MATLAB search path.

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* For questions, please contact Bastien Queste (b.queste at uea ac uk)