# %%
##Will Snyder
##Cambridge Psychiatry Brain Mapping Unit, NIMH Section on Developmental Neurogenomics
##March 2023
##will.snyder989@gmail.com

##The following script is a wrapper for computing sulcal phenotypes
##that are the basis for sulcal phenotype network (SPN) analysis
##All computations are performed on the sulcal graph and sulcal skeleton files. See
##brainvisa_helper_functions.py for how sulcal phenotypes are derived from these files
##See https://brainvisa.info/web/download.html for info on how to download and install BrainVISA

#import helper functions to compute phenotypes, as well as BrainVISA's image handling libraries
from brainvisa_helper_functions import *
from soma import aims
import pandas as pd
import shutil
import errno

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno in (errno.ENOTDIR, errno.EINVAL):
            shutil.copy(src, dst)
        else: raise

def extract_sulcal_phenotypes(fs_subject_name, output_dir, relabel_file, sulci_analyzed):
    print("Beginning sulcal phenotype measurement and SPN computation...")
    #open up list of sulci which BrainVISA reports morphometry on by default
    sulcus_file = open("/usr/local/spn_scripts/sulcus_list.txt")
    file_contents = sulcus_file.read()
    sulcus_list = file_contents.splitlines()
    #open up list of sulci to relabel, for which the first string
    #in each row is the applied label to all other sulci in that row
    label_retranslations = open(relabel_file)
    label_contents = label_retranslations.read()
    relabel_coding = label_contents.splitlines()

    #open up your sulcal graph and relabel it according to your relabel scheme
    morphologist_database = output_dir + '/morphologist_spn_database'
    left_graph_path =  morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/L" + fs_subject_name + "_default_session_auto.arg"
    right_graph_path = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/R" + fs_subject_name + "_default_session_auto.arg"
    left_data_path =  morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/L" + fs_subject_name + "_default_session_auto.data"
    right_data_path = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/R" + fs_subject_name + "_default_session_auto.data"
    left_graph_path_new =  morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/L" + fs_subject_name + "_default_session_auto_relabeled.arg"
    right_graph_path_new = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/R" + fs_subject_name + "_default_session_auto_relabeled.arg"
    left_data_path_new =  morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/L" + fs_subject_name + "_default_session_auto_relabeled.data"
    right_data_path_new = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/R" + fs_subject_name + "_default_session_auto_relabeled.data"
    
    #make a copy of the sulcal graph. that way, if you want to rerun with a different sulcal labeling scheme (or "parcellation"), you don't overwrite original brainvisa automated sulcal labeling
    copyanything(left_graph_path, left_graph_path_new)
    copyanything(right_graph_path, right_graph_path_new)
    copyanything(left_data_path, left_data_path_new)
    copyanything(right_data_path, right_data_path_new)
    left_sulcal_graph = aims.read(left_graph_path_new)
    right_sulcal_graph = aims.read(right_graph_path_new)
    print("relabeling sulcal graph according to your custom labeling scheme...")
    relabel_graph(left_sulcal_graph, relabel_coding, 'left',left_graph_path_new)
    relabel_graph(right_sulcal_graph, relabel_coding, 'right',right_graph_path_new)
    #CHANGE RELABELING SO YOU READ/WRITE TO A NEW GRAPH LOCATED AT A NEW PATH
    
    #read in the relabeled data
    left_relabeled_graph = aims.read(left_graph_path_new)
    right_relabeled_graph = aims.read(right_graph_path_new)

    #retrieve hull junction image (with sulcal exteriors labeled for each sulcal label),
    #sulcal skeleton image (which gives info on full sulcal medial wall) for processing in subsquent functions
    left_sulcal_skeleton_path = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/segmentation/Lskeleton_" + fs_subject_name + ".nii.gz"
    right_sulcal_skeleton_path = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/segmentation/Rskeleton_" + fs_subject_name + ".nii.gz"
    left_sulcal_skeleton = aims.read(left_sulcal_skeleton_path)
    right_sulcal_skeleton = aims.read(right_sulcal_skeleton_path)
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
    df.to_csv(morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/" + fs_subject_name + "_sulcal_phenotypes.csv")
    supp_table = pd.read_csv("/usr/local/spn_scripts/Snyder_Neuron_Table_S1.csv")
    group_efi = np.array(supp_table.iloc[:,2])
    
    sulcal_phenotypes = pd.read_csv(morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/" + fs_subject_name + "_sulcal_phenotypes.csv")
    compute_spn_and_efi(sulcal_phenotypes, sulci_analyzed, morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/" + fs_subject_name, group_efi, english_labels=True)
    
    print("Done! Sulcal phenotypes, SPN, end subject EFI each written as csv files")
