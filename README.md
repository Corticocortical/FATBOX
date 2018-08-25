# FATBOX
**The FMRI Analysis Toolbox**

The FATBOX (FMRI Analysis Toolbox) is a MATLAB and bash toolbox designed to simplify the analysis of fMRI data. It simplifies the following processes in particular: quickly creating ROIs in your volumetric scans, extracting data from ROIs (images or lists), integrating SPM, Freesurfer, and marsbar in your workflow, and interfacing with some multivariate tools (in particular RSA).

The functions (and scripts and examples) fall into a number of categories: 

**fMRI quality control**

These functions and scripts help you to assess the quality of your functional images, for instance by quantifying motion or the amplitude spectrum of your data. Note that some of these functions can build on top of each other by refining using each other's output as input (see respective documentation).

**fMRI temporal filtering**

Recent advances in fMRI protocol design have made it possible to scan at blisteringly fast acquisition speeds. If your acquisition speed is very fast (think < 1s), and your task is very slow (think long trials/mini-blocks at > ~ 12s), it may become useful to temporally filter your functional scans. Specifically, you may want to remove signal that is *faster than your task rate*. The assumption (!) here is that everything faster than your task rate is also faster than whatever process you are interested in, i.e. some form of noise.

In slow experiments with fast acquisition speeds (such as a multiband sequence), temporal filtering of the functional images can lead to respectable gains in sensitivity. However, it is best if you first test this using either an independent dataset or on a contrast that is orthogonal to your experimental manipulation. 

**General**

This is a collection of very generic functions that do nothing special. They are just convenient little helpers.
  
**Image handling**

These functions are mostly relevant to ROI and related types of analysis. These functions allow you to quickly create, inspect, and extract information from volumetric masks. There are also example scripts that show how marsbar and Freesurfer can be integrated in your workflow. For instance, this gives you the ability to create surface labels using Freesurfers sophisticated tools, and turn those surface labels into volumes that can conveniently be read out. Keep in mind that the example scripts are created for the computing environment at the Donders Institute.

**SPM control and wrappers**

These are functions that make it a little bit more convenient and faster to command SPM through the command line, such as to reslice, transform, segment, or apply contrasts to images. Note that many of these functions rely on the SPM batch system. When SPM updates, there might be compatibility issues. These functions should work with SPM8 and/or SPM12.

**Ext**

There are no files in this folder because these are functions and scripts from others. If a function in the FATBOX requires such functions, you can save them here.

**RSA**

This is a collection of function that simply representation similarity analysis. These functions are best to be used in conjunction with the RSA toolbox (see https://doi.org/10.1371/journal.pcbi.1003553). For that reason - and because RSA analysis is necessarily tied to the design of an experiment - some of the functions here are a little bit more idiosyncratic and might require some tailoring to suit your needs.
