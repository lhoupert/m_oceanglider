#!/bin/sh
## 
## Copyright (c) 2013 University of Washington.  All rights reserved.
## 
## This file contains proprietary information and remains the 
## unpublished property of the University of Washington. Use, disclosure,
## or reproduction is prohibited.
##
piloting_tarball="Piloting.tar"
linux_options="--owner 0 --group 0" 
# On Mac OS X, prevent writing resource fork files
COPYFILE_DISABLE=true
export COPYFILE_DISABLE
COPY_EXTENDED_ATTRIBUTES_DISABLE=true
export COPY_EXTENDED_ATTRIBUTES_DISABLE
macosx_options="--exclude='._*' --exclude='.svn' --exclude='.DS_Store' --exclude='*.bak' --exclude='*~'"
options=$macosx_options

# Start w/ documentation on dive process operation and coding observations
#tar $options -c -f $piloting_tarball README_piloting -C../docs Seaglider_Quality_Control_Manual.html
#2013.10.15, change.yohman@kongsberg.com, docs folder does not exist
tar $options -c -f $piloting_tarball README_piloting Seaglider_Quality_Control_Manual.html

# Diveplot graphs
tar $options -r -f $piloting_tarball diveplot.m diveplot_graphs.m
# Trim
tar $options -r -f $piloting_tarball trim.m
# Regress VBD, pitch, and roll support
tar $options -r -f $piloting_tarball regress_vbd.m w_misfit_func.m w_rms_func.m pitch_roll_control_bias.m c_line.m flightvec2.m
# Review results from piloting
tar $options -r -f $piloting_tarball review.m ts_diagram.m ts_diagram_dc.m waterfall_tsv.m

tar $options -r -f $piloting_tarball glide_slope.m flightvec.m find_stalled.m compute_oxygen_saturation.m cml_dens.m
# Bin a deployment
tar $options -r -f $piloting_tarball make_bin.m bin_edit.m bin_edit_dc.m bin_edit_dir.m bin_compare.m
# Read result files and unpack the results (ncload is for documentation purposes)
tar $options -r -f $piloting_tarball get_dive_data.m ncload.m unpack_data.m save_dive_results.m compute_dive_climb_indices.m compute_sensor_sg_indices.m interp1d.m unix_to_datenum.m datenum_to_unix.m
# Read and write log and eng files
tar $options -r -f $piloting_tarball parsecfg.m read_eng.m read_header.m read_log2.m read_gpctd.m
# Setup default constants
tar $options -r -f $piloting_tarball setup_constants.m clear_sg_config_constants.m sg_config_constants.m 
# Support QC operations (ala QC.py)
tar $options -r -f $piloting_tarball assert_qc.m bad_qc.m inherit_qc.m interpolate_data.m manual_qc.m qc_declarations.m

# Add miscellenanous support files
tar $options -r -f $piloting_tarball make_piloting_tarball.sh available_runs.m ask_which_runs.m available_profiles.m succinct_elts.m ctr1stdiffderiv.m  header.m plot_bounds.m underscore.m spice.m avail_mem.m

# billr: add some files that seem to have been overlooked
tar $options -r -f $piloting_tarball get_nc_data.m nc_bin_edit.m nc_bin_edit_dc.m pitch_regress.m read_sg_calib_constants.m regress_pitch.m regress_roll.m regress_vbd_dives.m reprocess.m rms.m sg_calib_constants.m sg_validate.m

# Seawater routines
tar $options -r -f $piloting_tarball sw_adtg.m sw_alpha.m sw_aonb.m sw_beta.m sw_bfrq.m sw_c3515.m sw_cndr.m sw_copy.m sw_cp.m sw_dens.m sw_dens0.m
tar $options -r -f $piloting_tarball sw_dist.m sw_dpth.m sw_f.m sw_fp.m sw_g.m sw_gpan.m sw_gvel.m sw_info.m sw_new.m sw_pden.m sw_pres.m
tar $options -r -f $piloting_tarball sw_ptmp.m sw_salds.m sw_salrp.m sw_salrt.m sw_sals.m sw_salt.m sw_seck.m sw_smow.m sw_svan.m sw_svel.m sw_temp.m
# compress it
gzip -f $piloting_tarball

# per Fritz:
# rsync -arv README_piloting Piloting.tar.gz glider@archimedes.ocean.washington.edu:.
