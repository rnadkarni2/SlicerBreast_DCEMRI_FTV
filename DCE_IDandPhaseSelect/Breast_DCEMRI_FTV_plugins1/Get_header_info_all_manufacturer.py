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
#script that contains all manufacturer specific
#functions for identifying the DCE series folders in an MR exam.


import pydicom
import dicom
import os
import re
import numpy as np
from statistics import mode
import hashlib

#Class FolderInfo is class that is used to construct objects with all of the header info for a particular folder in an exam.
#one class object per DICOM folder within exam
class FolderInfo:
    def __init__(self,folderpath):
        #switched from sitk method to os.listdir method because it is much faster
        #edit 5/29/2020: added verification that we're reading .dcm files because some numbered folders might not have DICOMs in them
        files = [f for f in os.listdir(folderpath) if f.endswith('.dcm')]
        FILES = [f for f in os.listdir(folderpath) if f.endswith('.DCM')]
        files_noext = [f for f in os.listdir(folderpath) if f.isdigit()] #edit 1/26/2021: In some folders, there is a series of DICOM images
                                                                         #with no .dcm or .DCM extension, but all the files have a number as the filename.

        if(len(files)>0 or len(FILES)>0 or len(files_noext) > 0):

            if len(files)>0:
                dcm_files = sorted(files)
                img1path = os.path.join(folderpath,dcm_files[0])
                imgendpath = os.path.join(folderpath,dcm_files[len(dcm_files)-1])

            if len(FILES)>0:
                dcm_files = sorted(FILES)
                img1path = os.path.join(folderpath,dcm_files[0])
                imgendpath = os.path.join(folderpath,dcm_files[len(dcm_files)-1])

            if len(files_noext)>0:
                dcm_files = sorted(files_noext)
                img1path = os.path.join(folderpath,dcm_files[0])
                imgendpath = os.path.join(folderpath,dcm_files[len(dcm_files)-1])

            print("Reading header from file:")
            print(img1path)

            try:
                header = pydicom.dcmread(img1path,stop_before_pixels = True)
                headerend = pydicom.dcmread(imgendpath,stop_before_pixels = True)
            except:
                try:
                    header = pydicom.dcmread(img1path,stop_before_pixels = True,force=True)
                    headerend = pydicom.dcmread(imgendpath,stop_before_pixels = True,force=True)
                except:
                    header = dicom.read_file(img1path)
                    headerend = dicom.read_file(imgendpath)

            #fill header info
            self.studydate = header.StudyDate

            #Edit 11/24/2020: Make totslices accurate for gunzipped cases
            if('gunzipped' in folderpath):
                allpath = folderpath.replace('\\gunzipped','')
            else:
                allpath = folderpath

            allfiles = os.listdir(allpath)
            self.totslices = len(allfiles) #Edit 11/2/2020: total # of slices in the folder, across all phases
            print("TOTAL SLICES")
            print(self.totslices)

            #Edit 10/30/2020: Include institution name field in this structure
            try:
                self.institution = header.InstitutionName
            except:
                try:
                    self.institution = header.ClinicalTrialSiteName
                except:
                    self.institution = 'Site Unknown'

            try:
                self.manufacturer = header.Manufacturer
            except:
                self.manufacturer = '' #Some images (MIPs?) do not have manufacturer field

            try:
                self.contenttime = round(float(header.ContentTime)) #Edit 1/20/2021: See if rounding to nearest whole number helps with error for 34306 v20
            except:
                self.contenttime = '000000'

            try:
                self.serdesc = header.SeriesDescription
            except:
                self.serdesc = ''


            #Temporal resolution in seconds based on manufacturer
            if ('GE' in self.manufacturer):
                try:
                    self.tempressec = header.TemporalResolution/1000
                except:
                    self.tempressec = 0

            if ('SIEMENS' in self.manufacturer or 'Siemens' in self.manufacturer):
                try:
                    self.tempressec = (header[0x19,0x100b].value)/1000
                except:
                    self.tempressec = 0

            #Acquisition duration for Philips
            if('Philips' in self.manufacturer or 'PHILIPS' in self.manufacturer):
                try:
                    self.acqdur = header.AcquisitionDuration
                #Edit 1/29/2021: try header[0x2005,0x1033].value when header.AcquisitionDuration doesn't work
                except:
                    try:
                        self.acqdur = header[0x2005,0x1033].value
                    #Edit 1/29/2021: Just set self.acqdur to 0 if neither of these work
                    except:
                        self.acqdur = 0

            try:
                self.lenimtype = len(header.ImageType) #Edit 7/15/2020: number of elements in image type array
                self.imtype0 = header.ImageType[0]
                self.imtype1 = header.ImageType[1]
                self.imtype2 = header.ImageType[2]

            except:
                self.lenimtype = 0
                self.imtype0 = ''
                self.imtype1 = ''
                self.imtype2 = ''

            try:
                self.slcorient = header[0x2001,0x100b].value
            except:
                self.slcorient = ''

            try:
                self.numtemp = header.NumberOfTemporalPositions
            except:
                self.numtemp = 0

            #may want to make this manufacturer dependent rather than multiple try-except
            try: #Siemens version
                self.pulseseqname = header.SequenceName
            except:
                try: #GE version
                    self.pulseseqname = header[0x19,0x109c].value

                    #edit 5/29/2020: convert from bytes to string if necessary
                    if isinstance(self.pulseseqname,bytes):
                        self.pulseseqname = self.pulseseqname.decode()

                except:
                    self.pulseseqname = ''

            try:
                self.matrix = header.AcquisitionMatrix
            except:
                self.matrix = [0,0,0,0]


            #Edit 6/1/2020: UCSF exam gave me bytes type error, although it doesn't have [21,104f] field.
            #For now, add exception to set numvol to 0
            try:
                self.im_in_acq = header[0x21,0x104f].value #using locations in acquisition field value
                #Edit 6/1/2020: Number of slices in an image. In some cases, folder contains entire DCE series, so # images in folder is some multiple of im_in_acq

                self.numvol = int(len(dcm_files)/(self.im_in_acq)) #Number of volumes is (# files in folder)/(# im_in_acq)
            except:
                #assume 1 volume in folder
                self.im_in_acq = len(dcm_files)
                self.numvol = 1

            try:
                self.pixsize0 = header.PixelSpacing[0]
                self.pixsize1 = header.PixelSpacing[1]
            except:
                self.pixsize0 = ''
                self.pixsize0 = ''

            try:
                self.orient = header.ImageOrientationPatient
            except:
                self.orient = ''

            try:
                self.tr = header.RepetitionTime
                self.trend = headerend.RepetitionTime
            except:
                self.tr = ''
                self.trend = ''

            #try filling in other group 18 DICOM fields
            try:
                self.numav = header.NumberOfAverages
            except:
                self.numav = ''

            try:
                self.pctsamp = header.PercentSampling
            except:
                self.pctsamp = ''

            try:
                self.SAR = header.SAR
            except:
                self.SAR = ''

            try:
                self.frame_ref = header.FrameOfReferenceUID
            except:
                self.frame_ref = ''

            try:
                self.fa = header.FlipAngle
                self.faend = headerend.FlipAngle
            except:
                self.fa = ''
                self.faend = ''

            #GE only
            try:
                self.fatsat = header[0x19,0x10a4].value
            except:
                self.fatsat = ''

            #Edit 11/23/2020: Transmit gain is used for GE folder identification sometimes
            try:
                self.transmitgain = round(header[0x19,0x10f9].value)
            except:
                self.transmitgain = ''

            try:
                self.transmitgainv2 = round(header[0x19,0x1094].value)
            except:
                self.transmitgainv2 = ''


            #Convert transmitgain from bytes if necessary
            if(isinstance(self.transmitgain,bytes)):
                self.transmitgain = int(self.transmitgain.decode())
            else:
                if(self.transmitgain != ''):
                    self.transmitgain = int(self.transmitgain)


            #Edit 2/23/21: Adding Image Position (Patient) to this structure
            #so that it can be added to matrix, im_in_acq, ... check
            try:
                self.imgpospat = header[0x20,0x32].value
            except:
                self.imgpospat = [0,0,0]

            #Edit 5/25/2021: Use private tags to find Manufacturer when
            #this information is not found in the default Manufacturer field.
            if('GE' not in self.manufacturer and 'PHILIPS' not in self.manufacturer and 'Philips' not in self.manufacturer and 'SIEMENS' not in self.manufacturer and 'Siemens' not in self.manufacturer):
                try:
                    privtag1 = header[0x19,0x10].value
                except:
                    privtag1 = ''

                try:
                    privtag2 = header[0x29,0x10].value
                except:
                    privtag2 = ''

                try:
                    privtag3 = header[0x2001,0x10].value
                except:
                    privtag3 = ''

                try:
                    privtag4 = header[0x2005,0x10].value
                except:
                    privtag4 = ''

                if('GE' in privtag1 or 'GE' in privtag2 or 'GE' in privtag3 or 'GE' in privtag4):
                    self.manufacturer = 'GE'

                if('SIEMENS' in privtag1 or 'SIEMENS' in privtag2 or 'SIEMENS' in privtag3 or 'SIEMENS' in privtag4):
                    self.manufacturer = 'SIEMENS'

                if('Siemens' in privtag1 or 'Siemens' in privtag2 or 'Siemens' in privtag3 or 'Siemens' in privtag4):
                    self.manufacturer = 'Siemens'

                if('PHILIPS' in privtag1 or 'PHILIPS' in privtag2 or 'PHILIPS' in privtag3 or 'PHILIPS' in privtag4):
                    self.manufacturer = 'PHILIPS'

                if('Philips' in privtag1 or 'Philips' in privtag2 or 'Philips' in privtag3 or 'Philips' in privtag4):
                    self.manufacturer = 'Philips'
                    print("manufacturer is Philips")


            #for GE, don't know what DICOM fields imlps and filter_mode are pulling info from
        #If folder doesn't have DICOM's, add warning of this as series description


def fillExamFolderInfoStructures(exampath):
    #return list of folders in exampath
    folders = [directory for directory in os.listdir(exampath) if os.path.isdir(os.path.join(exampath,directory))]

    #Edit 5/29/2020: Instead of using folder name as criterion, check if the folder actually has
    #DICOMs in and add to img_folders only if it does have DICOMs.
    img_folders = []
    for i in range(len(folders)):
        curr_path = os.path.join(exampath,folders[i])
        curr_files = [f for f in os.listdir(curr_path) if f.endswith('.dcm')]
        curr_FILES = [f for f in os.listdir(curr_path) if f.endswith('.DCM')]
        curr_files_noext = [f for f in os.listdir(curr_path) if f.isdigit()] #edit 1/26/2021: In some folders, there is a series of DICOM images
                                                                         #with no .dcm or .DCM extension, but all the files have a number as the filename.

        #Edit 7/6/2020: Program code to not count SlicerReports folder as an image folder
        #Edit 10/28/2020: Program code to only look at folders that have numerical names
        #Edit 7/6/2021: To make compatible with TCIA Duke exams, must consider folders that start
        #with a number as image folders too.
        if( (len(curr_files)>0 or len(curr_FILES)>0 or len(curr_files_noext)>0) and ('SlicerReport' not in folders[i]) and (folders[i].isdecimal() or folders[i][0].isdecimal()) ):
            if(folders[i].isdecimal()):
                img_folders.append(int(folders[i]))
            else:
                if(folders[i][0].isdecimal()):
                    img_folders.append(folders[i])

    img_folders = sorted(img_folders) #Edit 2/16/2021: Sort the folders in increasing numerical order
                                      #because that is more user-friendly for manual DCE folder ID.
    print(img_folders)

    #using list of DICOM image folders, make list of structures with DICOM header info from each folder
    all_folders_info = []
    for j in range(len(img_folders)):
        curr_path = os.path.join(exampath,str(img_folders[j]))
        curr_folder_info = FolderInfo(curr_path)
        all_folders_info.append(curr_folder_info)


    return img_folders, all_folders_info




#---------- Philips ----------#

#5/28/2020: Since David confirmed that pre and all post-contrast are almost always in the same folder, and 601 is not a
#pre-contrast folder for ISPY ID 20274, it appears that this is all I need to identify the DCE folder for Philips.

def philips_folder_lbl(img_folders, all_folders_info):

    #this is how to return all folder #'s that correspond to original, non-derived image (all folder #'s that have mod 1 when divided by 1)
    def condition(x): return (x%100) == 1
    original_folders = [element for idx, element in enumerate(img_folders) if condition(element)]
    original_idx = [idx for idx, element in enumerate(img_folders) if condition(element)]

    print("original folders")
    print(original_folders)

    #this is how to loop through only the original, non-derived images to find which one is 'dyn' (which has all DCE or post-contrast DCE only)
    dce_folders = []
    dce_ind = []
    for l in range(len(original_idx)):
        finf = all_folders_info[original_idx[l]]

        if ( ('dyn' in finf.serdesc) or ('DYN' in finf.serdesc) ):
            dce_folders.append(original_folders[l])
            dce_ind.append(original_idx[l])

    #if not identified by 'dyn' or 'DYN' in series description, use number of temporal positions > 1 to identify DCE folder
    if (len(dce_folders) == 0):
        #loop through all original folders to find one with number of temporal positions > 1
        for l in range(len(original_idx)):
            finf = all_folders_info[original_idx[l]]

            if ( finf.numtemp > 1 ):
                dce_folders.append(original_folders[l])
                dce_ind.append(original_idx[l])

    print("dce folders")
    print(dce_folders)

    return all_folders_info, dce_folders, dce_ind


#---------- Siemens ----------#


def siemens_folder_lbl(img_folders, all_folders_info):


    #Try to identify DCE images by presence of substring 'fl3d' or 'fl3_we' in pulseseqname field
    dce_folders = []
    dce_ind = [] #indices of dce folders in all_folders array
    t1_folders = []
    t1_ind = []
    for k in range(len(all_folders_info)):
        currfinf = all_folders_info[k]

        #Edit 7/15/2020: If (0008,0008) Image Type field has more than 4 elements, series is derived and shouldn't be added to original DCE list
        #Edit 7/17/2020: Remove image type length restriction because it seems to be causing errors in folder identification
        #and it didn't help with removing MIPs from DCE list anyway
        #Edit 7/22/2020: Only count as DCE if first field in ImageType field is ORIGINAL
        if ( (('fl3d' in currfinf.pulseseqname) or ('fl3_we' in currfinf.pulseseqname)) and ('ORIGINAL' in currfinf.imtype0) ):
            dce_folders.append(img_folders[k])
            dce_ind.append(k)

    print("DCE folders initial guess based on pulse sequence info")
    print(dce_folders)

    #Remove images with SUB in series description from the DCE list and add them to subtraction list
    dce_folders_orig = dce_folders
    dce_ind_orig = dce_ind
    dce_folders = []
    dce_ind = []

    #Edit 7/15/2020: Make sure that SUB and MIP are both not in series description
    dce_submip_folders = []
    dce_submip_ind = []
    for s in range(len(dce_ind_orig)):
        finf = all_folders_info[dce_ind_orig[s]]
        #Edit 2/2/21: Add 'sub' to this check due to UMinn 19984 v30
        if( ('SUB' not in finf.serdesc) and ('sub' not in finf.serdesc) and ('MIP' not in finf.serdesc) ):
            dce_folders.append(dce_folders_orig[s])
            dce_ind.append(dce_ind_orig[s])
        else:
            dce_submip_folders.append(dce_folders_orig[s])
            dce_submip_ind.append(dce_ind_orig[s])

    print("dce folders after removal of subtractions and MIPs")
    print(dce_folders)

    #Edit 2/5/21: Make this matrix, im_in_acq ... check more complicated, with max of 5 groups instead of max of 2.

    folders1 = []
    idx1 = []

    folders2 = []
    idx2 = []

    folders3 = []
    idx3 = []

    folders4 = []
    idx4 = []

    folders5 = []
    idx5 = []


    #Function to check for match in matrix, im_in_acq ... match
    def paramCheckFunc(foldersn, idxn, fold_inf, a, match_found):
        foldersn_0info = all_folders_info[idxn[0]]

        #2/24/2021: Added image position patient to this check for Siemens too.
        if (fold_inf.matrix == foldersn_0info.matrix and fold_inf.im_in_acq == foldersn_0info.im_in_acq and fold_inf.pixsize0 == foldersn_0info.pixsize0 and fold_inf.pixsize1 == foldersn_0info.pixsize1 and fold_inf.orient == foldersn_0info.orient and fold_inf.fatsat == foldersn_0info.fatsat and fold_inf.imgpospat[0] == foldersn_0info.imgpospat[0] and fold_inf.imgpospat[1] == foldersn_0info.imgpospat[1] and fold_inf.imgpospat[2] == foldersn_0info.imgpospat[2]):
            foldersn.append(dce_folders[a])
            idxn.append(dce_ind[a])
            match_found = 1

        return foldersn, idxn, match_found

    match_found = 0 #variable to check when to add to new group
    num_group = 0 #variable to keep track of number of groups divided by common matrix, im_in_acq, ...
    print("checking for match in matrix im_in_acq, pixsizes, orient, fatsat")
    for a in range(len(dce_ind)):
        match_found = 0 #reset to 0 at start of each iteration
        fold_inf = all_folders_info[dce_ind[a]]

        if a == 0:
            #If first folder in DCE list, just add its info to folders1 and idx1
            folders1.append(dce_folders[a])
            idx1.append(dce_ind[a])
            num_group = 1
            match_found = 1
        else:
            #First, check for parameter match with folders1
            folders1, idx1, match_found = paramCheckFunc(folders1, idx1, fold_inf, a, match_found)

            #Then, move to folders2
            if(len(folders2) > 0):
                if(match_found == 0):
                    #if folders 2 is not empty and match not found yet, check for parameter match
                    folders2, idx2, match_found = paramCheckFunc(folders2, idx2, fold_inf, a, match_found)
            else:
                #if match not found yet and folders2 and idx2 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders2.append(dce_folders[a])
                    idx2.append(dce_ind[a])
                    num_group = 2
                    match_found = 1


            #Then, move to folders3
            if(len(folders3) > 0 ):
                if(match_found == 0):
                    #if folders 3 is not empty and match not found yet, check for parameter match
                    folders3, idx3, match_found = paramCheckFunc(folders3, idx3, fold_inf, a, match_found)
            else:
                #if match not found yet and folders3 and idx3 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders3.append(dce_folders[a])
                    idx3.append(dce_ind[a])
                    num_group = 3
                    match_found = 1

            #Then, move to folders4
            if(len(folders4) > 0 ):
                if(match_found == 0):
                    #if folders 4 is not empty and match not found yet, check for parameter match
                    folders4, idx4, match_found = paramCheckFunc(folders4, idx4, fold_inf, a, match_found)
            else:
                #if match not found yet and folders4 and idx4 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders4.append(dce_folders[a])
                    idx4.append(dce_ind[a])
                    num_group = 4
                    match_found = 1


            #Then, move to folders5
            if(len(folders5) > 0 ):
                if(match_found == 0):
                    #if folders 5 is not empty and match not found yet, check for parameter match
                    folders5, idx5, match_found = paramCheckFunc(folders5, idx5, fold_inf, a, match_found)
            else:
                #if match not found yet and folders5 and idx5 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders5.append(dce_folders[a])
                    idx5.append(dce_ind[a])
                    num_group = 5
                    match_found = 1

    print("different groups based on matrix, im_in_acq, pixsizes, orient, fatsat")
    print(folders1)
    print(folders2)
    print(folders3)
    print(folders4)
    print(folders5)

    #Use loop to assign DCE as the group with max # of totslices
    totslicesmax = 0 #variable to keep track of which group has most totslices
    for zzzz in range(num_group):
        totslices_zzzz = 0 #number of totslices for group #zzzz

        if(zzzz == 0):
            f = folders1
            idx = idx1

        if(zzzz == 1):
            f = folders2
            idx = idx2

        if(zzzz == 2):
            f = folders3
            idx = idx3

        if(zzzz == 3):
            f = folders4
            idx = idx4

        if(zzzz == 4):
            f = folders5
            idx = idx5

        #In this loop, sum up all the totslices for group #zzzz
        for yyyy in range(len(idx)):
            infof = all_folders_info[idx[yyyy]]
            totslices_zzzz = totslices_zzzz + infof.totslices

        #If we have achieved a new max, assing dce to this and update totslicesmax
        if(totslices_zzzz > totslicesmax):
            totslicesmax = totslices_zzzz
            dce_folders = f
            dce_ind = idx


    print("dce folders after matrix, im_in_acq ... check")
    print(dce_folders)

    #Filter out some values from DCE based on TR
    #Assume the mode TR value is the DCE one
    trvals = []
    for aa in range(len(dce_ind)):
        fold_inf = all_folders_info[dce_ind[aa]]
        trvals.append(fold_inf.tr)

    #If you can get a trmode, use it to narrow down which folders are dce
    try:
        trmode = mode(trvals)
        dce_folders_trmode = []
        dce_ind_trmode = []
        for bb in range(len(dce_ind)):
            fold_inf = all_folders_info[dce_ind[bb]]

            #if tr value matches mode,it is dce. otherwise, it is t1
            if( fold_inf.tr == trmode ):
                dce_folders_trmode.append(dce_folders[bb])
                dce_ind_trmode.append(dce_ind[bb])
            else:
                t1_folders.append(dce_folders[bb])
                t1_ind.append(dce_ind[bb])

        dce_folders = dce_folders_trmode
        dce_ind = dce_ind_trmode
        print("DCE folders after T_R check")
        print(dce_folders)
    except:
        print("Unable to use T_R for dce_folder identification")

    #Edit 1/27/2021: Must make sure that all folders labeled as DCE have the same number of slices
    #Edit 2/8/2021: Restrict this check to > 4 dce_folders
    if(len(dce_folders)>4):
        alltotslices = []
        for aaa in range(len(dce_ind)):
            curr_finf = all_folders_info[dce_ind[aaa]]
            alltotslices.append(curr_finf.totslices)

        #Edit 2/9/2021: If only 2 DCE folders with different # of totslices, just choose higher totslices value
        #This relies on the assumption that Siemens will never have 2 DCE folders, but will have either 1 or >= 3.
        if( len(dce_ind) == 2 and alltotslices[0] != alltotslices[1] ):
            totslices_mode = max(alltotslices)
            totslices_modecount = 1
        else:
            totslices_mode = mode(alltotslices)
            totslices_modecount = alltotslices.count(totslices_mode) #Edit 4/20/21: return # of times mode occurs

        #If totslices in a folder is not equal to mode value of totslices, delete that folder from list of DCE folders
        bbb = 0
        while (bbb < len(dce_ind)):
            curr_finf = all_folders_info[dce_ind[bbb]]

            if(curr_finf.totslices != totslices_mode):
                #4/20/21: If number of slices for current folder exceeds total # of slices
                #associated with totslices_mode, current folder is the only one that is DCE.
                if(curr_finf.totslices > totslices_mode*totslices_modecount):
                    dce_folders = [dce_folders[bbb]]
                    dce_ind = [dce_ind[bbb]]
                    break #exit while loop
                #Otherwise, remove current folder from DCE list
                else:
                    del(dce_folders[bbb])
                    del(dce_ind[bbb])
            else:
                bbb = bbb + 1

    #Once you have the final list of DCE folders, sort them in order of increasing Content Time
    ctimes = []
    for z in range(len(dce_ind)):
        finf = all_folders_info[dce_ind[z]]
        ctimes.append(finf.contenttime)
    ctimes = np.array(ctimes)
    ctind = ctimes.argsort()
    ctimes = ctimes[ctind]

    dce_folders = np.array(dce_folders)
    dce_folders = dce_folders[ctind]
    dce_folders = dce_folders.tolist()

    dce_ind = np.array(dce_ind)
    dce_ind = dce_ind[ctind]
    dce_ind = dce_ind.tolist()

    print("Content times sorted")
    print(ctimes)

    print("dce folders sorted by content time")
    print(dce_folders)

    print("t1 folders")
    print(t1_folders)

    return all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind, dce_submip_folders, dce_submip_ind


#---------- GE ----------#
#Comment 6/1/2020: This UCSF exam doesn't have Locations in Acquisition or Images in Acquisition DICOM fields. What to do?
def ge_folder_lbl(img_folders, all_folders_info):

    #Try to identify DCE images by presence of substring 'efgre3d' or '3dfgre' in pulseseqname field
    dce_folders = []
    dce_ind = [] #indices of dce folders in all_folders array
    for k in range(len(all_folders_info)):
        currfinf = all_folders_info[k]
        #Edit 7//15/2020: Pulse sequence has ispy in its name for some Columbia exams
        #Edit 7/22/2020: Must have ORIGINAL in 1st element of ImageType field to be considered DCE
        if ( (('efgre3d' in currfinf.pulseseqname) or ('3dfgre' in currfinf.pulseseqname) or ('ispy' in currfinf.pulseseqname)) and ('ORIGINAL' in currfinf.imtype0) ):
            dce_folders.append(img_folders[k])
            dce_ind.append(k)

    print("DCE folders initial guess based on pulse sequence alone")
    print(dce_folders)
    print("dce ind")
    print(dce_ind)
    print("image folders all")
    print(img_folders)

    #Try to remove CAD (processed) images from dce folders list
    non_cad_folders = []
    non_cad_ind = []
    print("Running section to remove CAD images from DCE folder list")
    for c in range(len(dce_ind)):
        finf = all_folders_info[dce_ind[c]]

        #Edit 1/19/21: Adding PRE-SCAN in series description as a factor to eliminate from DCE folder list
        if( ('SUB' not in finf.serdesc) and ('SAG' not in finf.serdesc) and ('COR' not in finf.serdesc) and ('MIP' not in finf.serdesc) and ('PRE-SCAN' not in finf.serdesc) and ('SCOUT' not in finf.serdesc) ):
            non_cad_folders.append(dce_folders[c])
            non_cad_ind.append(dce_ind[c])

    if( len(non_cad_folders) > 0):
        dce_folders = non_cad_folders
        dce_ind = non_cad_ind

    print("dce folders after removal of CAD images")
    print(dce_folders)
    print("dce ind")
    print(dce_ind)
    print("image folders all")
    print(img_folders)

    #Edit 7/15/2020: If dce_folders array is empty after scanning through pulse sequence descriptions, use series descriptions
    #to find out which folders are DCE
    if( len(dce_folders) == 0):

        #Edit 1/22/2021: Added 'or's to these if statements so that it works for UCSD 36061 v20 as well.

        #Edit 7/20/2020: First, just check for 'Ph' to find all post-contrast series. Then, you know that the part of the series description
        #that comes after 'Ph#/' is the pre-contrast series description, so you should use this to find the pre-contrast folder
        for ki in range(len(all_folders_info)):
            currfinf = all_folders_info[ki]

            #5/20/2021: Check for other identifiers to make sure it is axial, non-subtraction post-contrast image
            if('Ph' in currfinf.serdesc or 'post' in currfinf.serdesc and ('SUB' not in currfinf.serdesc) and ('SAG' not in currfinf.serdesc) and ('COR' not in currfinf.serdesc) and ('MIP' not in currfinf.serdesc) and ('PRE-SCAN' not in currfinf.serdesc) and ('SCOUT' not in currfinf.serdesc) ):
                dce_folders.append(img_folders[ki])
                dce_ind.append(ki)

        #4/26/21: Only attempt to look for pre-contrast image if post-contrast images
        #were found based on series description
        if len(dce_folders) > 0:
            ex_post_folder_info = all_folders_info[dce_ind[0]] #take folder info for one of the post-contrast folders
            ex_post_serdesc = ex_post_folder_info.serdesc #take series description of this post-contrast folder
            pre_serdesc = ex_post_serdesc[4:] #use post-contrast series description to obtain pre-contrast series description
            print("Pre-contrast should have this series description:")
            print(pre_serdesc)
            #loop through folders again to find pre-contrast one
            for li in range(len(all_folders_info)):
                currfinf = all_folders_info[li]
                #find pre-contrast folder and append it to your dce list
                if( (pre_serdesc in currfinf.serdesc and 'Ph' not in currfinf.serdesc) or 'pre' in currfinf.serdesc ): #Need 'Ph' not in to avoid having 2 copies of each post-contrast series
                    dce_folders.append(img_folders[li])
                    dce_ind.append(li)
        print("dce folders after checking series description")
        print(dce_folders)

    #Try to remove images with (#/#/#)-(#/#/#) in series description, because these are subtraction images
    non_sub_folders = []
    non_sub_ind = []
    print("Check for (#/#/#)-(#/#/#) series description corresponding to subtraction images")

    for s in range(len(dce_ind)):
        finf = all_folders_info[dce_ind[s]]

        numbers = re.findall('[0-9]+',finf.serdesc) #find all numbers in series description of current folder DICOMs
        slashinds = [i for i, letter in enumerate(finf.serdesc) if letter == '/'] #find all slashes in current series description
        if ( len(numbers) < 6 or len(slashinds) < 4): #at least 6 numbers and at least 4 / in series description likely corresponds to (#/#/#)-(#/#/#)
            non_sub_folders.append(dce_folders[s])
            non_sub_ind.append(dce_ind[s])

    if(len(non_sub_folders)>0):
        dce_folders = non_sub_folders
        dce_ind = non_sub_ind

    print("DCE folders after removal of (#/#/#)-(#/#/#) subtraction images")
    print(dce_folders)

    #Try to remove image folders in which the image doesn't have the same flip angle and repetition time in first and last slice
    dce_fatrmatch_folders = []
    dce_fatrmatch_ind = []
    print("Checking for mismatch between flip angle and repetition time in first and last slice")

    for m in range(len(dce_ind)):

        finf = all_folders_info[dce_ind[m]]
        if (finf.tr == finf.trend and finf.fa == finf.faend):
            dce_fatrmatch_folders.append(dce_folders[m])
            dce_fatrmatch_ind.append(dce_ind[m])

    if len(dce_fatrmatch_folders) > 0:
        dce_folders = dce_fatrmatch_folders
        dce_ind = dce_fatrmatch_ind

    print("dce folders after check for flip angle and repetition time match")
    print(dce_folders)

    #Out of the folders that have pulse sequence name corresponding to DCE, find the ones with 'ORIG' in series
    #description because these are the original, non-processed ones

    orig_folders = []
    orig_ind = []
    print("checking for ORIG in series description")
    for o in range(len(dce_ind)):
        finf = all_folders_info[dce_ind[o]]

        if ('ORIG' in finf.serdesc):
            orig_folders.append(dce_folders[o])
            orig_ind.append(dce_ind[o])

    #If folders were found with ORIG in DICOM series descriptions, these are the dce folders we want, so set dce variables accordingly
    if (len(orig_folders)>0):
        dce_folders = orig_folders
        dce_ind = orig_ind

    print("dce folders after check for ORIG in series description")
    print(dce_folders)

    #Edit 2/3/2021: Delete folders with < 32 slices from DCE folder list
    zzz = 0
    while zzz < len(dce_ind):
        finfzzz = all_folders_info[dce_ind[zzz]]

        if(finfzzz.totslices < 32):
            del(dce_folders[zzz])
            del(dce_ind[zzz])
        else:
            zzz = zzz+1



    #From DCE folder list, detect any erroneous addition of T1 based on
    #difference in matrix, im_in_acq, pixsizes, orient, gesat, imlps,
    #flip_first, trs_first, filter_mode

    #Edit 6/1/2020: idl code says allow for small variations in flip angle or TR. Not sure what values to use
    #for small variation, so just excluding these parameters from difference checks for now
    #2nd Edit 6/1/2020: You are supposed to assume that there is only one combination of parameters for max # of volumes,
    #and that the folders that have this parameter combo and max # of volumes are (for now) DCE folders
    #numvolmax is being incorporated because of the GE exams that have entire DCE series in one folder

    #Edit 6/1/2020: Revert back to non numvolmax dependent method of choosing dce and t1 folders for now
    ###Construct two empty list. In the end, one will contain only dce folder #'s
    ###and the other will contain only T1 folder #'s

    #Edit 2/5/21: Make this matrix, im_in_acq ... check more complicated, with max of 5 groups instead of max of 2.

    folders1 = []
    idx1 = []

    folders2 = []
    idx2 = []

    folders3 = []
    idx3 = []

    folders4 = []
    idx4 = []

    folders5 = []
    idx5 = []


    #Function to check for match in matrix, im_in_acq ... match
    def paramCheckFunc(foldersn, idxn, fold_inf, a, match_found):
        foldersn_0info = all_folders_info[idxn[0]]

        print("transmit gains")
        print(fold_inf.transmitgain)
        print(fold_inf.transmitgainv2)
        print("")
        print(foldersn_0info.transmitgain)
        print(foldersn_0info.transmitgainv2)


        #Edit 2/8/2021: Added transmitgain and transmitgainv2 to this check
        #Edit 2/23/2021: Added image position patient to this check
        if (fold_inf.matrix == foldersn_0info.matrix and fold_inf.im_in_acq == foldersn_0info.im_in_acq and fold_inf.pixsize0 == foldersn_0info.pixsize0 and fold_inf.pixsize1 == foldersn_0info.pixsize1 and fold_inf.orient == foldersn_0info.orient and fold_inf.fatsat == foldersn_0info.fatsat and  fold_inf.transmitgain == foldersn_0info.transmitgain and fold_inf.transmitgainv2 == foldersn_0info.transmitgainv2 and fold_inf.imgpospat[0] == foldersn_0info.imgpospat[0] and fold_inf.imgpospat[1] == foldersn_0info.imgpospat[1] and fold_inf.imgpospat[2] == foldersn_0info.imgpospat[2] ):
            foldersn.append(dce_folders[a])
            idxn.append(dce_ind[a])
            match_found = 1

        return foldersn, idxn, match_found

    match_found = 0 #variable to check when to add to new group
    num_group = 0 #variable to keep track of number of groups divided by common matrix, im_in_acq, ...
    print("checking for match in matrix im_in_acq, pixsizes, orient, fatsat")
    for a in range(len(dce_ind)):
        match_found = 0 #reset to 0 at start of each iteration
        fold_inf = all_folders_info[dce_ind[a]]

        if a == 0:
            #If first folder in DCE list, just add its info to folders1 and idx1
            folders1.append(dce_folders[a])
            idx1.append(dce_ind[a])
            num_group = 1
            match_found = 1
        else:
            #First, check for parameter match with folders1
            folders1, idx1, match_found = paramCheckFunc(folders1, idx1, fold_inf, a, match_found)

            #Then, move to folders2
            if(len(folders2) > 0):
                if(match_found == 0):
                    #if folders 2 is not empty and match not found yet, check for parameter match
                    folders2, idx2, match_found = paramCheckFunc(folders2, idx2, fold_inf, a, match_found)
            else:
                #if match not found yet and folders2 and idx2 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders2.append(dce_folders[a])
                    idx2.append(dce_ind[a])
                    num_group = 2
                    match_found = 1


            #Then, move to folders3
            if(len(folders3) > 0 ):
                if(match_found == 0):
                    #if folders 3 is not empty and match not found yet, check for parameter match
                    folders3, idx3, match_found = paramCheckFunc(folders3, idx3, fold_inf, a, match_found)
            else:
                #if match not found yet and folders3 and idx3 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders3.append(dce_folders[a])
                    idx3.append(dce_ind[a])
                    num_group = 3
                    match_found = 1

            #Then, move to folders4
            if(len(folders4) > 0 ):
                if(match_found == 0):
                    #if folders 4 is not empty and match not found yet, check for parameter match
                    folders4, idx4, match_found = paramCheckFunc(folders4, idx4, fold_inf, a, match_found)
            else:
                #if match not found yet and folders4 and idx4 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders4.append(dce_folders[a])
                    idx4.append(dce_ind[a])
                    num_group = 4
                    match_found = 1


            #Then, move to folders5
            if(len(folders5) > 0 ):
                if(match_found == 0):
                    #if folders 5 is not empty and match not found yet, check for parameter match
                    folders5, idx5, match_found = paramCheckFunc(folders5, idx5, fold_inf, a, match_found)
            else:
                #if match not found yet and folders5 and idx5 are empty, populate these with the current dce folder
                if(match_found == 0):
                    folders5.append(dce_folders[a])
                    idx5.append(dce_ind[a])
                    num_group = 5
                    match_found = 1

    print("different groups based on matrix, im_in_acq, pixsizes, orient, fatsat")
    print(folders1)
    print(folders2)
    print(folders3)
    print(folders4)
    print(folders5)

    #Use loop to assign DCE as the group with max # of totslices
    totslicesmax = 0 #variable to keep track of which group has most totslices
    for zzzz in range(num_group):
        totslices_zzzz = 0 #number of totslices for group #zzzz

        if(zzzz == 0):
            f = folders1
            idx = idx1

        if(zzzz == 1):
            f = folders2
            idx = idx2

        if(zzzz == 2):
            f = folders3
            idx = idx3

        if(zzzz == 3):
            f = folders4
            idx = idx4

        if(zzzz == 4):
            f = folders5
            idx = idx5

        #In this loop, sum up all the totslices for group #zzzz
        for yyyy in range(len(idx)):
            infof = all_folders_info[idx[yyyy]]
            totslices_zzzz = totslices_zzzz + infof.totslices

        #If we have achieved a new max, assing dce to this and update totslicesmax
        if(totslices_zzzz > totslicesmax):
            totslicesmax = totslices_zzzz
            dce_folders = f
            dce_ind = idx

    print("dce folders after matrix, im_in_acq ... check")
    print(dce_folders)

    #Edit 11/3/2020: Folder ident fix that is specific to UCSD ID 04967 v30
    #If there are 2 DCE folders and one has 'Post' in the name and the other still has more image slices, the other is DCE.
    if(len(dce_ind) == 2):
        dce1info = all_folders_info[dce_ind[0]]
        dce2info = all_folders_info[dce_ind[1]]

        if('Post' in dce1info.serdesc and dce2info.totslices > dce1info.totslices):
            dce_ind = [dce_ind[1]]
            dce_folders = [dce_folders[1]]

        if('Post' in dce2info.serdesc and dce1info.totslices > dce2info.totslices):
            dce_ind = [dce_ind[0]]
            dce_folders = [dce_folders[0]]

    t1_folders = []
    t1_ind = []

    if(len(dce_ind) > 1):
        #Filter out some values from DCE based on TR
        #Assume the mode TR value is the DCE one
        trvals = []
        totslicevals = []
        for aa in range(len(dce_ind)):
            fold_inf = all_folders_info[dce_ind[aa]]
            trvals.append(fold_inf.tr)
            totslicevals.append(fold_inf.totslices)

        trvals_unique = np.unique(trvals)

        if(len(trvals)>len(trvals_unique)):
            #use trmode if you can
            try:
                trmode = mode(trvals)
                dce_folders_trmode = []
                dce_ind_trmode = []
                for bb in range(len(dce_ind)):
                    fold_inf = all_folders_info[dce_ind[bb]]

                    #if tr value matches mode,it is dce. otherwise, it is t1
                    if( fold_inf.tr == trmode ):
                        dce_folders_trmode.append(dce_folders[bb])
                        dce_ind_trmode.append(dce_ind[bb])
                    else:
                        t1_folders.append(dce_folders[bb])
                        t1_ind.append(dce_ind[bb])

                dce_folders = dce_folders_trmode
                dce_ind = dce_ind_trmode
            except:
                print("Each folder currently labeled as DCE doesn't have unique T_R, but we're unable to use T_R for dce folder identification")
        else:
            #Edit 1/21/2021: If each folder in current DCE list has unique T_R, only the one with greatest # of totslices is DCE
            ind_maxslc = np.argmax(totslicevals)
            dce_folders = [dce_folders[ind_maxslc]]
            dce_ind = [dce_ind[ind_maxslc]]

    print("t1 folders after T_R check")
    print(t1_folders)

    #From pair of folders, see if they have different group 18 values
    #if they do, only the one with more phases in the folder is DCE.
    print("Check for difference in group 18 values")
    if(len(dce_ind) == 2):
        dce1info = all_folders_info[dce_ind[0]]
        dce2info = all_folders_info[dce_ind[1]]

        if(round(dce1info.tr,3) != round(dce2info.tr,3) or round(dce1info.numav,3) != round(dce2info.numav,3) or round(dce1info.pctsamp,3) != round(dce2info.pctsamp,3)  or round(dce1info.SAR,3) != round(dce2info.SAR,3) ):
            print("The 2 folders are not both DCE. Choosing one to be the DCE folder.")
            #if the 2 folders are not matching in all these parameters, pick the one that has greater # phases. #phase = totslices/im_in_acq
            if( (dce1info.totslices)/(dce1info.im_in_acq) > dce2info.totslices/(dce2info.im_in_acq)):
                prntstr1 = "Folder #" + str(dce_folders[0]) + " is DCE."
                print(prntstr1)
                dce_ind = [dce_ind[0]]
                dce_folders = [dce_folders[0]]
            else:
                prntstr2 = "Folder #" + str(dce_folders[1]) + " is DCE."
                print(prntstr2)
                dce_ind = [dce_ind[1]]
                dce_folders = [dce_folders[1]]

    print("DCE folders after group 18 check")
    print(dce_folders)

    #For UCSD, need to exclude series with 'TEST' in description from dce-- add this to t1 folders list
    #Edit 5/29/2020: Changing this so that it doesn't always add last folder to t1 list
    b = 0
    while b < (len(dce_ind)):
        fold_inf = all_folders_info[int(dce_ind[b])]

        if ('TEST' in fold_inf.serdesc):
            t1_folders.append(dce_folders[b])
            t1_ind.append(dce_ind[b])
            del dce_folders[b]
            del dce_ind[b]
        else:
            b = b+1

    print("DCE folders after removing TEST images")
    print(dce_folders)


    #Edit 2/8/2021: Keep this transmit gain check even though that check has been added to the matrix, im_in_acq, ... check.

    #From pair of folders, see if they have different transmit gain
    #if they do, only the one with more phases in the folder is DCE.
    print("Check for transmit gain difference")
    if(len(dce_ind) == 2):
        dce1info = all_folders_info[dce_ind[0]]
        dce2info = all_folders_info[dce_ind[1]]

        #Edit 1/15/2021: Try different way of using both header fields with transmit gain
        if(dce1info.transmitgain != dce2info.transmitgain or dce1info.transmitgainv2 != dce2info.transmitgainv2):
            print("The 2 folders are not both DCE. Choosing one to be the DCE folder.")
            #if the 2 folders are not matching in all these parameters, pick the one that has greater # phases. #phase = totslices/im_in_acq
            if( ((dce1info.totslices)/(dce1info.im_in_acq)) >  ((dce2info.totslices)/(dce2info.im_in_acq)) ):
                prntstr1 = "Folder #" + str(dce_folders[0]) + " is DCE."
                print(prntstr1)
                dce_ind = [dce_ind[0]]
                dce_folders = [dce_folders[0]]
            else:
                prntstr2 = "Folder #" + str(dce_folders[1]) + " is DCE."
                print(prntstr2)
                dce_ind = [dce_ind[1]]
                dce_folders = [dce_folders[1]]
    print("DCE folders after transmit gain check")
    print(dce_folders)


    #Edit 1/28/2021: Add this same # of folders check to ge_folder_lbl. Only do it for > 4 DCE folders
    #because GE can have 1, 2, 5, 6, or 7 DCE folders

    #Edit 1/27/2021: Must make sure that all folders labeled as DCE have the same number of slices
    if(len(dce_folders) > 4):
        alltotslices = []
        for aaa in range(len(dce_ind)):
            curr_finf = all_folders_info[dce_ind[aaa]]
            alltotslices.append(curr_finf.totslices)
        totslices_mode = mode(alltotslices)
        totslices_modecount = alltotslices.count(totslices_mode) #Edit 4/20/21: return # of times mode occurs

        #If totslices in a folder is not equal to mode value of totslices, delete that folder from list of DCE folders
        bbb = 0
        while (bbb < len(dce_ind)):
            curr_finf = all_folders_info[dce_ind[bbb]]

            if(curr_finf.totslices != totslices_mode):
                #4/20/21: If totslices for current folder is > total # of slices
                #associated with mode value of totslices, only current folder is DCE.
                if(curr_finf.totslices > totslices_mode*totslices_modecount):
                    dce_folders = [dce_folders[bbb]]
                    dce_ind = [dce_ind[bbb]]
                    break #exit loop
                #Otherwise, delete current folder from DCE list
                else:
                    del(dce_folders[bbb])
                    del(dce_ind[bbb])
            else:
                bbb = bbb + 1

    #Once you have the final list of DCE folders, sort them in order of increasing Content Time
    ctimes = []
    for z in range(len(dce_ind)):
        finf = all_folders_info[dce_ind[z]]
        ctimes.append(finf.contenttime)
    ctimes = np.array(ctimes)
    ctind = ctimes.argsort()
    ctimes = ctimes[ctind]

    dce_folders = np.array(dce_folders)
    dce_folders = dce_folders[ctind]
    dce_folders = np.unique(dce_folders)
    dce_folders = dce_folders.tolist()

    dce_ind = np.array(dce_ind)
    dce_ind = dce_ind[ctind]
    dce_ind = np.unique(dce_ind)
    dce_ind = dce_ind.tolist()

    print("Content times sorted")
    print(ctimes)

    print("dce folders sorted by content time")
    print(dce_folders)

    print("t1 folders")
    print(t1_folders)

    return all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind
