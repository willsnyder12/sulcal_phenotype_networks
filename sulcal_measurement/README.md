The files sulcal_phenotype_wrapper.py and brainvisa_helper_functions.py comprise the sulcal phenotype extraction pipeline.
The output of running sulcal_phenotype_wrapper.py is sulcal_phenotypes.csv, giving a value to each of the 5 sulcal phenotypes
for all sulci measureable. The wrapper script is run within a directory containing the sulcal graph ({left/right}_sulcal_graph.arg and {left/right}_sulcal_graph.data files/directories) and sulcal skeleton ({left/right}_sulcal_skeleton.nii) outputs from BrainVISA. See
brainvisa_helper_functions.py for how sulcal phenotypes are derived from these files.

Wrapper code is run by copying all SPN / sulcal phenotype processing files into the subject's sulcal graph directory, then using a singularity command like
`singularity run -B ${SULCAL_GRAPH_DIRECTORY}:${SULCAL_GRAPH_DIRECTORY} /path/to/your/brainvisa/singularity/file/brainvisa-5.0.4.sif python sulcal_phenotype_wrapper.py`
See https://brainvisa.info/web/download.html for info on how to download and install BrainVISA. Note that trimesh must be installed by entering into BrainVISA's bash distribution (e.g., `bv bash`) and using `pip install trimesh`

Prior to publication, a full singularity pipeline will be available that automatically outputs sulcal graph and skeleton files with the correct file names for this pipeline, as well as automatically extracting sulcal phenotypes thereafter, allowing a single line of code to take you from FreeSurfer outputs to sulcal_phenotypes.csv.

Credit for the function for circumscribed circle area (smallestenclosingcircle.py) is https://www.nayuki.io/page/smallest-enclosing-circle
Credit for the function for fractal fimension calculation (FractalDimension.py) is https://github.com/ChatzigeorgiouGroup/FractalDimension
