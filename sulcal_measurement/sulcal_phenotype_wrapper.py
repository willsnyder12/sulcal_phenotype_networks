# %%
##Will Snyder
##Cambridge Psychiatry Brain Mapping Unit, NIMH Section on Developmental Neurogenomics
##March 2023
##will.snyder989@gmail.com

##The following script is a wrapper for computing sulcal phenotypes
##that are the basis for sulcal phenotype network (SPN) analysis
##Code must be run using BrainVISA's python distribution, within a directory containig sulcal graphs (*.arg files)
##(and their corresponding *.data folders) output by BrainVISA's Morphologist. Also trimesh will also need
##to be installed within BrainVISA's python distribution with `pip install trimesh`
##All computations are performed on the sulcal graph and sulcal skeleton files. See
##brainvisa_helper_functions.py for how sulcal phenotypes are derived from these files

##Wrapper code is run by copying all SPN / sulcal phenotype processing files into the subject's sulcal graph
##directory, then using a singularity command like
##`singularity run -B ${SULCAL_GRAPH_DIRECTORY}:${SULCAL_GRAPH_DIRECTORY} /path/to/your/brainvisa/singularity/file/brainvisa-5.0.4.sif python sulcal_phenotype_wrapper.py`
##See https://brainvisa.info/web/download.html for info on how to download and install BrainVISA

#import helper functions to compute phenotypes, as well as BrainVISA's image handling libraries
from brainvisa_helper_functions  import *
from soma import aims
import pandas as pd

#open up list of sulci which BrainVISA reports morphometry on by default
sulcus_file = open("sulcus_list.txt")
file_contents = sulcus_file.read()
sulcus_list = file_contents.splitlines()
#open up list of sulci to relabel, for which the first string
#in each row is the applied label to all other sulci in that row
label_retranslations = open("relabel_list.txt")
label_contents = label_retranslations.read()
relabel_coding = label_contents.splitlines()

#open up your sulcal graph and relabel it according to your relabel scheme
left_graph_path = "left_sulcal_graph.arg"
right_graph_path = "right_sulcal_graph.arg"
left_sulcal_graph = aims.read(left_graph_path)
right_sulcal_graph = aims.read(right_graph_path)
print("relabeling sulcal graph according to your custom labeling scheme...")
relabel_graph(left_sulcal_graph, relabel_coding, 'left',left_graph_path)
relabel_graph(right_sulcal_graph, relabel_coding, 'right',right_graph_path)

#read in the relabeled data
left_relabeled_graph = aims.read(left_graph_path)
right_relabeled_graph = aims.read(right_graph_path)

#retrieve hull junction image (with sulcal exteriors labeled for each sulcal label),
#sulcal skeleton image (which gives info on full sulcal medial wall) for processing in subsquent functions
left_sulcal_skeleton = aims.read("left_sulcal_skeleton.nii")
right_sulcal_skeleton = aims.read("right_sulcal_skeleton.nii")
print("creating sulcus hull junction images for feature extraction...")
left_hj = write_hull_junctions(left_relabeled_graph,left_sulcal_skeleton,sulcus_list, 'left')
right_hj = write_hull_junctions(right_relabeled_graph,right_sulcal_skeleton,sulcus_list, 'right')


##compute sulcal phenotypes
##RADIAL MEASURES
#get depth profiles to extract depth / "radial" measures
print("computing geodesic depth profiles...")
geodesic_depth_left = compute_depth_profiles(left_relabeled_graph, left_sulcal_skeleton, sulcus_list)
geodesic_depth_right = compute_depth_profiles(right_relabeled_graph, right_sulcal_skeleton, sulcus_list)

##TANGENTIAL MEASURES
#compute longest branch length
print("computing maximum length of each sulcus' connected components...")
left_longest_branch = compute_largest_cc_length(left_relabeled_graph,sulcus_list)
right_longest_branch = compute_largest_cc_length(right_relabeled_graph,sulcus_list)
#compute branch span
print("computing sulcus branch span...")
branch_span_left = compute_branch_span(left_hj,left_sulcal_skeleton,sulcus_list)
branch_span_right = compute_branch_span(right_hj,left_sulcal_skeleton,sulcus_list)
#compute fd
print("computing fractal dimensoon...")
fd_left = compute_fractal_dimension(left_hj,sulcus_list)
fd_right = compute_fractal_dimension(right_hj,sulcus_list)

#Combine across hemispheres, and compute median and median absolute deviation of depth profiles
print("combining hemisphere features...")
average_depth, depth_variability, longest_branch, branch_span, fd = combine_hemi_features(sulcus_list,
                                                                                    geodesic_depth_left,
                                                                                    geodesic_depth_right,
                                                                                    left_longest_branch,
                                                                                    right_longest_branch,
                                                                                    branch_span_left,
                                                                                    branch_span_right,
                                                                                    fd_left,
                                                                                    fd_right)
df = pd.DataFrame({'sulcus_label': sulcus_list,'average_depth':average_depth,'depth_variability':depth_variability,
                   'longest_branch':longest_branch,'branch_span':branch_span, 'fd': fd})
print(df)
#write out result as a csv file using default BrainVISA sulcal labels, showing no values where you 
#defined a sulcal label to be merged with another sulcal label, and also showing NA if no data was available
#i.e. no branches of that sulcus were identified by BrainVISA Morphologist
df.to_csv("sulcal_phenotypes.csv")
print("Done! Sulcal phenotypes written as csv file")





