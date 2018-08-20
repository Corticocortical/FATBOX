#!/bin/tcsh/
#
# Example to use  Freesurfer through Torque.
#
# This file can be submitted using qsub.
# To do so, you must be in the UNIX T C Shell (type tcsh).
# You should then set the Freesurfer home directory, the subject 
# directory and source the Freesurfer setup, with the commands
# commented out below.
#
# USE THESE TO SET UP FREESURFER:
# setenv FREESURFER_HOME /opt/freesurfer/5.3
# setenv SUBJECTS_DIR /home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects
# source $FREESURFER_HOME/SetUpFreeSurfer.csh
#
# SUBMIT THE RECONSTRUCTION:
# qsub Reconstruction.csh -l walltime=10:00:00,mem=5gb
#
# With this file being a script, you can theoretically also submit all subjects
# at once, by making the subject ID (and potentially the reconstruction mode)
# arguments to pass to this script. The script would recognize these as $1, and
# $2, respectively. Then you might be able to just call this script 24 times, 
# although I have not verified that qsub actually works that way.
recon-all -subjid S40 -all
