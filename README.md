# Sulcal Phenotype Networks (SPNs)
Code for SPN analysis as described in "A bipolar taxonomy of adult human brain sulcal morphology related to timing of fetal sulcation and trans-sulcal gene expression gradients"

The folder sulcal_measurement contains python code for extracting the 5 sulcal phenotypes (average depth, depth variability, longest branch length, branch span, and fractal dimension)
using BrainVISA's python distribution.

The folder spn_analysis contains R code for the network analysis using these measures, highlighting main figures and results.

The folder sulcal_parcellations contains sulcal parcellations used in surface-based analyses.

Prior to publication, a full singularity pipeline will be available that automatically outputs sulcal graph and skeleton files with the correct file names for this pipeline, as well as automatically extracting sulcal phenotypes thereafter, allowing a single line of code to take you from FreeSurfer outputs to a csv file of sulcal phenotypes and a csv file for the SPN matrix.
