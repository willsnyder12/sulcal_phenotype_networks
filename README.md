# Sulcal Phenotype Networks (SPNs)
Code for SPN analysis as described in "A bimodal taxonomy of adult human brain sulcal morphology related to timing of fetal sulcation and trans-sulcal gene expression gradients"
Please find the containerized version of this pipeline at: 

dx.doi.org/10.6084/m9.figshare.25874425

This container takes as input a FreeSurfer subject directory and outputs BrainVISA sulcal morphometry + sulcal phenotypes + sulcal phenotype networks + sulcal complexity measurement (eigenfold index)

The folder spn_scripts contains code used in the container for extracting the 5 sulcal phenotypes (average depth, depth variability, longest branch length, branch span, and fractal dimension) in part using BrainVISA's python distribution.

The folder spn_analysis contains R code for the network analysis using these measures, highlighting main figures and results.

The folder sulcal_parcellations contains sulcal parcellations used in surface-based analyses.

The folder visualization contains additional scripts to visualize results in R.

If you run into any questions, feel free to email me @ will.snyder989@gmail.com !

<img width="600" alt="Screenshot 2024-08-22 at 2 47 20 PM" src="https://github.com/user-attachments/assets/8eb73ce5-d1bf-495d-a3ae-ffb5248ab2ce">





