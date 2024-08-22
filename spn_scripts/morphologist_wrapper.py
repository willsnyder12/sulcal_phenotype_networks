from brainvisa import axon
from brainvisa.axon import processes as axon_processes
import os
import sys
import shutil
import errno

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno in (errno.ENOTDIR, errno.EINVAL):
            shutil.copy(src, dst)
        else: raise

def run_morphologist_from_freesurfer(fs_subject_name, fs_dir, output_dir): 
    #initialize brainvisa processes
    axon.initializeProcesses()
    context = axon_processes.defaultContext()

    #if first time running morphologist pipeline, set up databases for importing / managing inputs and outputs
    morphologist_database = output_dir + '/morphologist_spn_database'
    freesurfer_database = output_dir + '/freesurfer_spn_database'
    if not os.path.exists(morphologist_database):
        print("Setting up freesurfer and brainvisa databases in specified output directory...")
        db_proc = axon_processes.getProcessInstance('Create Database')
        db_proc.database_directory = morphologist_database
        db_proc.ontology = 'brainvisa-3.2.0'
        context.runProcess(db_proc)
        db_proc.database_directory = freesurfer_database
        db_proc.ontology = 'freesurfer'
        context.runProcess(db_proc)

    #copy input freesurfer data into the freesurfer spn database, if you haven't already
    fs_subject_database_dir = freesurfer_database + "/" + fs_subject_name
    if not os.path.exists(fs_subject_database_dir):
        copyanything(fs_dir + "/" + fs_subject_name, fs_subject_database_dir)

    #import t1 to brainvisa morphologist database
    bv_import_proc = axon_processes.getProcessInstance('Import T1 MRI')
    bv_import_proc.input = fs_subject_database_dir + "/mri/orig.mgz"
    bv_import_proc.output_database = morphologist_database
    context.runProcess(bv_import_proc)

    #CHANGE SO THAT MULTITHREADING IS TURNED OFF
    fs_import_proc = axon_processes.getProcessInstance('Import FreeSurfer grey/white segmentation to Morphologist')
    fs_import_proc.T1_orig = fs_subject_database_dir + "/mri/orig.mgz"
    fs_import_proc.ribbon_image = fs_subject_database_dir + "/mri/ribbon.mgz"
    fs_import_proc.Talairach_Auto = fs_subject_database_dir + "/mri/transforms/talairach.auto.xfm"
    fs_import_proc.T1_output = morphologist_database + "/subjects/" + fs_subject_name + "/t1mri/default_acquisition/" + fs_subject_name + ".nii.gz"
    fs_import_proc.allow_multithreading = False
    context.runProcess(fs_import_proc)
