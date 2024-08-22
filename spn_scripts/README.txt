The files sulcal_phenotype_wrapper.py and brainvisa_helper_functions.py comprise the main code for the sulcal phenotype extraction pipeline. The output from running these in the container is a csv file, giving a value to each of the 5 sulcal phenotypes for all 40 sulci measured. Additionally the 40x40 within-subject sulcal phenotype network (SPN) is computed.

ComputeSPN is the wrapper used in the container, which based on inputs sends info to spn_pipeline.py, which then runs BrainVISA Morphologist and then extracts sulcal phenotypes.

Credit for the function for circumscribed circle area (smallestenclosingcircle.py) is https://www.nayuki.io/page/smallest-enclosing-circle Credit for the function for fractal fimension calculation (FractalDimension.py) is https://github.com/ChatzigeorgiouGroup/FractalDimension
The trimesh library is also used in brainvisa_helper_functions.py, with trimesh source code found at https://github.com/mikedh/trimesh
