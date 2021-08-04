#---------------------------------------------------------------------------------
#Copyright (c) 2021
#By Rohan Nadkarni and the Regents of the University of California.
#All rights reserved.

#This code was developed by the UCSF Breast Imaging Research Group,
#and it is part of the extension Breast_DCEMRI_FTV,
#which can be installed from the Extension Manager in the 3D Slicer application.

#If you intend to make derivative works of any code from the
#Breast_DCEMRI_FTV extension, please inform the
#UCSF Breast Imaging Research Group
#(https://radiology.ucsf.edu/research/labs/breast-imaging-research).

#This notice must be attached to all copies, partial copies,
#revisions, or derivations of this code.
#---------------------------------------------------------------------------------


#Created by Rohan Nadkarni
#function to check if an exam's DICOMs need
#to be gunzipped or not

import os


def checkForGzippedDicoms(exampath):
    #return list of folders in exampath
    folders = [directory for directory in os.listdir(exampath) if os.path.isdir(os.path.join(exampath,directory))]

    #by default, gzipped binary variable is set to 0
    gzipped = 0

    #loop through all folders to see if you find a gzipped DICOM
    #If you find one, exam is considered gzipped

    #Edit 12/16/2020: Change the criterion to be finding > 1 gzipped DICOMs
    for i in range(len(folders)):
        curr_path = os.path.join(exampath,folders[i])

        gzip_files_lc = [f for f in os.listdir(curr_path) if f.endswith('.dcm.gz')] #lowercase dcm gzipped
        gzip_files_uc = [f for f in os.listdir(curr_path) if f.endswith('.DCM.gz')] #lowercase dcm gzipped
        #Edit 6/30/21: make sure gzipped .dmi files are not counted in this list
        gzip_files_noext = [f for f in os.listdir(curr_path) if(f.endswith('.gz') and '.dmi' not in f)] #Edit 1/26/21: For DICOMs that don't have .dcm or .DCM extension

        #If current folder has either of these gzipped DICOM types,
        #considered exam to be gzipped and exit loop
        #6/30/21: Folder with gzipped exam must have number as name
        #This prevents gunzipping prompted by folders containing
        #gzipped .dmi files, for example.
        if(folders[i].isdigit() and (len(gzip_files_lc) > 1) or (len(gzip_files_uc) > 1) or (len(gzip_files_noext) > 1) ):
            gzipped = 1
            break

    return gzipped
