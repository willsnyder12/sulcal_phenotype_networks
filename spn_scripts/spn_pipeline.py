from morphologist_wrapper import *
from sulcal_phenotype_wrapper import *
import sys


n_args = len(sys.argv) - 1
if n_args == 3:
    fs_subject_name = sys.argv[1]
    fs_dir = sys.argv[2]
    output_dir = sys.argv[3]
    run_morphologist_from_freesurfer(fs_subject_name, fs_dir, output_dir)
    #give default values for sulcal relabeling file / file for sulci to be analyzed
    extract_sulcal_phenotypes(fs_subject_name, output_dir, "/usr/local/spn_scripts/relabel_list.txt", "/usr/local/spn_scripts/sulci_for_analysis.txt")
if n_args == 5:
    fs_subject_name = sys.argv[1]
    fs_dir = sys.argv[2]
    output_dir = sys.argv[3]
    relabel_file = sys.argv[4]
    sulci_analyzed = sys.argv[5]
    run_morphologist_from_freesurfer(fs_subject_name, fs_dir, output_dir)
    extract_sulcal_phenotypes(fs_subject_name, output_dir, relabel_file, sulci_analyzed)
else:
    raise RuntimeError("incorrect number of input arguments")