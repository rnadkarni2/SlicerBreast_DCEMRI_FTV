#---------------------------------------------------------------------------------
#Copyright (c) 2021
#By Rohan Nadkarni and the Regents of the University of California.
#All rights reserved.

#This code was developed by the UCSF Breast Imaging Research Group,
#and it is part of the extension Breast_DCEMRI_TumorSegment,
#which can be installed from the Extension Manager in the 3D Slicer application.

#If you intend to make derivative works of any code from the
#Breast_DCEMRI_TumorSegment extension, please inform the
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

##        reader = sitk.ImageSeriesReader()
##        dicom_names = reader.GetGDCMSeriesFileNames(folderpath)

##        img1path = dicom_names[0] #already includes full path
##        imglastpath = dicom_names[len(dicom_names)-1] #already includes full path

        #switched from sitk method to os.listdir method because it is much faster
        #edit 5/29/2020: added verification that we're reading .dcm files because some numbered folders might not have DICOMs in them
        files = [f for f in os.listdir(folderpath) if f.endswith('.dcm')]
        FILES = [f for f in os.listdir(folderpath) if f.endswith('.DCM')]
        files_noext = [f for f in os.listdir(folderpath) if f.isdigit()] #edit 1/26/2021: In some folders, there is a series of DICOM images
                                                                         #with no .dcm or .DCM extension, but all the files have a number as the filename.

    
        if(len(files)>0 or len(FILES)>0 or len(files_noext) > 0):


            if len(files)>0:
                dcm_files = sorted(files)
                img1path = folderpath + "\\" + dcm_files[0]
                imgendpath = folderpath + "\\" + dcm_files[len(dcm_files)-1]

            if len(FILES)>0:
                dcm_files = sorted(FILES)
                img1path = folderpath + "\\" + dcm_files[0]
                imgendpath = folderpath + "\\" + dcm_files[len(dcm_files)-1]

            if len(files_noext)>0:
                dcm_files = sorted(files_noext)
                img1path = folderpath + "\\" + dcm_files[0]
                imgendpath = folderpath + "\\" + dcm_files[len(dcm_files)-1]

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
##        else:
##            self.serdesc = 'NOT DICOM'
        


def fillExamFolderInfoStructures(exampath):
    #return list of folders in exampath
    folders = [directory for directory in os.listdir(exampath) if os.path.isdir(exampath + "\\" + directory)]

    #Edit 5/29/2020: Instead of using folder name as criterion, check if the folder actually has
    #DICOMs in and add to img_folders only if it does have DICOMs.
    img_folders = []
    for i in range(len(folders)):
        curr_path = exampath + "\\" + folders[i]
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
        curr_path = exampath + "\\" + str(img_folders[j])
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
    
    #print folder info structure for all image folders (original or otherwise)
##    for finf in all_folders_info:
##        print(finf.serdesc)
##        print(finf.imtype0)
##        print(finf.imtype1)
##        print(finf.slcorient)
##        print(finf.numtemp)
##        print("")
##    print(" ")
##    print(" ")
##    print(" ")

    return all_folders_info, dce_folders, dce_ind

##exampath = r"\\researchfiles.radiology.ucsf.edu\birp_perm2\ispy_2019\4460_OHSC\66493\20190129_v10\E20190129"
##print("---------- Philips exam identification ----------")
##img_folders, all_folders_info = fillExamFolderInfoStructures(exampath)
##all_folders_info, dce_folders, dce_ind = philips_folder_lbl(img_folders, all_folders_info)
##











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

    #print("dce subtraction or mip folders")
    #print(dce_submip_folders)



    ##print("dce folders")
    ##print(dce_folders)



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




##### OLD VERSION of matrix, im_in_acq, ... check. #####

    #From DCE folders list, detect any erroneous addition of 'T1' based on
    #difference in matrix, im_in_acq, pixel sizes, orient, frame_ref
    #Construct two empty list. In the end, one will contain only dce folder #'s
    #and the other will contain only T1 folder #'s
##    folders1 = []
##    idx1 = []
##
##    folders2 = []
##    idx2 = []
##    for a in range(len(dce_ind)):
##        fold_inf = all_folders_info[dce_ind[a]]
##        if a == 0:
##            matrix0 = fold_inf.matrix
##            im_in_acq0 = fold_inf.im_in_acq
##            pixsize00 = fold_inf.pixsize0
##            pixsize10 = fold_inf.pixsize1
##            frame_ref0 = fold_inf.frame_ref
##            orient0 = fold_inf.orient
##    
##        #This assumes there will only be 2 possible combos of these parameters, one for DCE and other for T1. I don't know if that is the case.
##        if (fold_inf.matrix == matrix0 and fold_inf.im_in_acq == im_in_acq0 and fold_inf.pixsize0 == pixsize00 and fold_inf.pixsize1 == pixsize10 and fold_inf.frame_ref == frame_ref0 and fold_inf.orient == orient0):
##            folders1.append(dce_folders[a])
##            idx1.append(dce_ind[a])
##        else:
##            folders2.append(dce_folders[a])
##            idx2.append(dce_ind[a])
##
##    #assume folders list with more elements is dce, and other one is T1
##    #Edit 2/5/2021: Add condition that len(folders1) must be <= 8 to be considered DCE.
##    if (len(folders1)>len(folders2)):
##
##        if(len(folders1) <= 8):
##            dce_folders = folders1
##            dce_ind = idx1
##    
##            t1_folders = folders2
##            t1_ind = idx2
##        else:
##            dce_folders = folders2
##            dce_ind = idx2
##    
##            t1_folders = folders1
##            t1_ind = idx1
##
##    #Edit 2/5/21: Instead of 'else', make this mirror the above 'if' structure
##    if (len(folders2)>len(folders1)):
##
##        if(len(folders2) <= 8):
##            dce_folders = folders2
##            dce_ind = idx2
##    
##            t1_folders = folders1
##            t1_ind = idx1
##        else:
##            dce_folders = folders1
##            dce_ind = idx1
##    
##            t1_folders = folders2
##            t1_ind = idx2


    

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


    #print folder info structure for all image folders (original or otherwise)
##    for finf in all_folders_info:
##        print(finf.serdesc)
##        print(finf.imtype0)
##        print(finf.imtype1)
##        print(finf.imtype2)
##        print(finf.slcorient)
##        print(finf.numtemp)
##        print(finf.pulseseqname)
##        print(finf.tr)
##        print("")
##
##    print(" ")
##    print(" ")
##    print(" ")

    return all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind, dce_submip_folders, dce_submip_ind

##exampath = r"\\researchfiles.radiology.ucsf.edu\birp_perm2\ispy_2019\961_Gtown\89396\20191002_v10\E9620302"
##print("---------- Siemens exam identification ----------")
##all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind, dce_sub_folders, dce_sub_ind = siemens_folder_lbl(exampath)












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

    #Edit 1/28/2021: While loop to remove folders that have T1 in series description
    #Edit 2/4/2021: Get rid of this section because the DCE folders 700-706 for SSWTI 18059 v30 have 'T1' in the series description
##    ccc = 0
##    while(ccc < len(dce_folders)):
##        finfccc = all_folders_info[dce_ind[ccc]]
##
##        if('T1' in finfccc.serdesc):
##            del dce_folders[ccc]
##            del dce_ind[ccc]
##        else:
##            ccc = ccc+1
            


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


    
    #Edit 11/16/2020: Series Description - based correction for UCSF ISPY ID 11749 v30
    
    #2/9/2021: I found out that 'T1' in Series Description doesn't mean the folder is not DCE,
    #so I'll comment out this section and hope that other edits will fix issues for exams that
    #used this section.
    
##    if(len(dce_folders) == 3):
##        f1 = all_folders_info[dce_ind[0]]
##        f2 = all_folders_info[dce_ind[1]]
##        f3 = all_folders_info[dce_ind[2]]
##
##        if('multiPhase' in f1.serdesc and 'T1' in f2.serdesc and 'T1' in f3.serdesc):
##            t1_folders = [img_folders[dce_ind[1]]]
##            t1_ind = [dce_ind[1]]
##            t1_folders.append(img_folders[dce_ind[2]])
##            t1_ind.append(dce_ind[2])
##
##            dce_folders = [img_folders[dce_ind[0]]]
##            dce_ind = [dce_ind[0]]
##            
##        if('multiPhase' in f2.serdesc and 'T1' in f1.serdesc and 'T1' in f3.serdesc):
##            t1_folders = [img_folders[dce_ind[0]]]
##            t1_ind = [dce_ind[0]]
##            t1_folders.append(img_folders[dce_ind[2]])
##            t1_ind.append(dce_ind[2])
##
##            dce_folders = [img_folders[dce_ind[1]]]
##            dce_ind = [dce_ind[1]]
##
##        if('multiPhase' in f3.serdesc and 'T1' in f1.serdesc and 'T1' in f2.serdesc):
##            t1_folders = [img_folders[dce_ind[0]]]
##            t1_ind = [dce_ind[0]]
##            t1_folders.append(img_folders[dce_ind[1]])
##            t1_ind.append(dce_ind[1])
##
##            dce_folders = [img_folders[dce_ind[2]]]
##            dce_ind = [dce_ind[2]]

    #Edit 11/16/2020: Series Description - based correction for UCSF ISPY ID 10936 v10
    #Edit 12/14/2020: Moved CAD removal to beginning and made the 'and's 'or's to fix folder ident
    #for MDAnd 12124

    #2/9/2021: I found out that 'T1' in Series Description doesn't mean the folder is not DCE,
    #so I'll comment out this section and hope that other edits will fix issues for exams that
    #used this section.

            
##    if(len(dce_folders) == 2):
##        print("Checking for multiPhase or T1 in series description")
##        
##        f1 = all_folders_info[dce_ind[0]]
##        f2 = all_folders_info[dce_ind[1]]
##
##        #Edit 1/19/2021: Edited these 'if' statements to also include 'not in'
##
##        if( ('multiPhase' in f1.serdesc and 'multiPhase' not in f2.serdesc) or ('T1' in f2.serdesc and 'T1' not in f1.serdesc) ):
##            print("folder 1 is DCE")
##            t1_folders = [img_folders[dce_ind[1]]]
##            t1_ind = [dce_ind[1]]
##
##            dce_folders = [img_folders[dce_ind[0]]]
##            dce_ind = [dce_ind[0]]
##            
##        if( ('multiPhase' in f2.serdesc and 'multiPhase' not in f1.serdesc) or ('T1' in f1.serdesc and 'T1'not in f2.serdesc) ):
##            print("folder 2 is DCE")
##            t1_folders = [img_folders[dce_ind[0]]]
##            t1_ind = [dce_ind[0]]
##
##            print(dce_ind[1])
##            print(img_folders)
##
##            dce_folders = [img_folders[dce_ind[1]]]
##            dce_ind = [dce_ind[1]]
##
##    print("DCE folders after multiphase and T1 check")
##    print(dce_folders)
##    print("dce ind")
##    print(dce_ind)
##    print("image folders all")
##    print(img_folders)

    

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
                
        

    ##print("Folders with DCE pulse sequence")
    ##print(dce_folders)




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


    ##print("Folders with DCE pulse sequence that are not CAD processed images")
    ##print(dce_folders)

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

    ##print("Folders with DCE pulse sequence that are not CAD processed and have gone through ORIG label test")
    ##print(dce_folders)
    

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

    #for numvolmax = 1 or 2, assume folders list with more elements is dce, and other one is T1, because 5-7 folders for DCE series
    #for numvolmax > 2, only one of the folders should have elements at all


    #Edit 6/1/2020: Revert back to non numvolmax dependent method of choosing dce and t1 folders for now

    #case where numvolmax > 2
    #this is for case where DCE phases are all in one folder
    ##if numvolmax > 2:
    ##    if (len(folders2) == 0):
    ##        dce_folders = folders1
    ##        dce_ind = idx1
    ##    else:
    ##        if (len(folders1) == 0):
    ##            dce_folders = folders2
    ##            dce_ind = idx2
    ##        else:
    ##            "Could not find DCE folders"
    #case where numvolmax is 2 or 1
    #since locations in acquisition can sometimes have 1/2 the number of slices even when
    #each DCE phase has its own folder, need to include value 2 as well.
    ##else:

    #Edit 6/1/2020: folder array with more elements is DCE, unless
    #the shorter folder array has PRE/POST in at least 1 series description,
    #indicating all DCE phases are in that folder
    #This PRE/POST method came from my own observations, not David's code


    #Edit 2/5/2021: Out of the 5 groups the one with the most slices is the DCE group
    

    ##

######################   OLD CODE for matrix, im_in_acq, ... CHECK   ################################

    
##    t1_folders = []
##    t1_ind = []
##    if len(folders1)>len(folders2):
##
##        exception = 0
##        for p in range(len(idx2)):
##            finf = all_folders_info[idx2[p]]
##            
##            finf_idx1 = all_folders_info[idx1[0]]
##            #Edit 1/28/2021: If > 1 folders in idx1, incorporate 2nd folder into totslices check too.
##            if(len(idx1) > 1):
##                finf_idx1_2 = all_folders_info[idx1[1]]
##            else:
##                finf_idx1_2.totslices = 0
##
##            #Edit 1/26/2021: Incorporate totslices as an exception. If group with only 1 folder has more slices per folder, this one folder is the DCE.
##            #Edit 1/29/2021: Incorporate # of folders in idx1 as an additional component of totslices check <-- CANCELING THIS EDIT BECAUSE IT CREATES ERRORS FOR OTHER EXAMS
##            if('PRE/POST' in finf.serdesc or 'PRE AND POST' in finf.serdesc or (finf.totslices > finf_idx1.totslices and finf.totslices > finf_idx1_2.totslices) ):
##                exception = 1
##
##        if(exception == 1):
##            dce_folders = folders2
##            dce_ind = idx2
##    
##            t1_folders = folders1
##            t1_ind = idx1
##        else:
##            dce_folders = folders1
##            dce_ind = idx1
##    
##            t1_folders = folders2
##            t1_ind = idx2
##        
##    else:
##        #Edit 11/24/2020: Case where you have 2 folders and have to use matrix/orient/fatsat check to separate other folder from multivolume DCE folder
##        if(len(folders1) == 1 and len(folders2) == 1):
##            f1information = all_folders_info[idx1[0]]
##            f2information = all_folders_info[idx2[0]]
##            
##            if(f1information.totslices/f1information.im_in_acq > 1):
##                print("folder 1 is multivolume DCE, folder 2 is not DCE")
##                dce_folders = folders1
##                dce_ind = idx1
##                
##                t1_folders = folders2
##                t1_ind = idx2
##                
##            if(f2information.totslices/f2information.im_in_acq > 1):
##                print("folder 2 is multivolume DCE, folder 1 is not DCE")
##                dce_folders = folders2
##                dce_ind = idx2
##
##                t1_folders = folders1
##                t1_ind = idx1
##            
##        else:
##            exception = 0
##            for p in range(len(idx1)):
##                finf = all_folders_info[idx1[p]]
##                finf_idx2 = all_folders_info[idx2[0]]
##
##                #Edit 1/28/2021: If > 1 folders in idx2, incorporate 2nd folder into totslices check too.
##                if(len(idx2) > 1):
##                    finf_idx2_2 = all_folders_info[idx2[1]]
##                else:
##                    finf_idx2_2.totslices = 0
##
##
##                #Edit 1/26/2021: Incorporate totslices as an exception. If group with only 1 folder has more slices per folder, this one folder is the DCE.
##                #Edit 1/29/2021: Incorporate # of folders in idx2 as an additional component of totslices check <-- CANCELING THIS EDIT BECAUSE IT CREATES ERRORS FOR OTHER EXAMS
##                if('PRE/POST' in finf.serdesc or 'PRE AND POST' in finf.serdesc or (finf.totslices > finf_idx2.totslices and finf.totslices > finf_idx2_2.totslices) ):
##                    exception = 1
##
##            if(exception == 1):
##                dce_folders = folders1
##                dce_ind = idx1
##    
##                t1_folders = folders2
##                t1_ind = idx2
##            else:
##                dce_folders = folders2
##                dce_ind = idx2
##    
##                t1_folders = folders1
##                t1_ind = idx1

##    print("t1 folders after matrix/orient/fatsat check")
##    print(t1_folders)
##
##    print("dce folders after check for match in im_in_acq, pixsizes, orient, fatsat")
##    print(dce_folders)


    
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
    
    ##print("t1folders after")
    ##print(t1_folders)

    #Edit 5/29/2020: Add T1 folder back to dce list if it has Ph in series description (eg Ph1, Ph2, etc)
    #Edit 6/3/2020: Comment out this section, because sometimes you have ORIG Ph1 and separate Ph1 folder, etc
    
##    t = 0
##    while ( t <= (len(t1_ind) -1) ):
##        finf = all_folders_info[int(t1_ind[t])]
##    
##        if('Ph' in finf.serdesc):
##            dce_folders.append(int(t1_folders[t]))
##            dce_ind.append(int(t1_ind[t]))
##
##            del t1_ind[t]
##            del t1_folders[t]
##
##        else:
##            t = t+1

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


    #print folder info structure for all image folders (original or otherwise)
##    for finf in all_folders_info:
##        print(finf.serdesc)
##        print(finf.imtype0)
##        print(finf.imtype1)
##        print(finf.imtype2)
##        print(finf.slcorient)
##        print(finf.numtemp)
##        print(finf.pulseseqname)
##        print(finf.tr)
##        print(finf.fa)
##        print(finf.fatsat)
##        print("")

    #Edit 10/30 2020: If UCSF and <= 3 DCE folders, the last DCE folder is the only true DCE folder.
##    folder1info = all_folders_info[0]
##    if('UCSF' in folder1info.institution and len(dce_folders) <= 3):
##        dce_ind = dce_ind[-1]
##        dce_ind = [dce_ind]
##        
##        dce_folders = dce_folders[-1]
##        dce_folders = [dce_folders]

    return all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind

#Call the GE folder label function for specific exam path
##exampath = r"\\researchfiles.radiology.ucsf.edu\birp_perm2\ispy_2019\4218_RMH\79922\20190828_v10\E32596"
##print("---------- GE exam identification ----------")
##all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind = ge_folder_lbl(exampath)




















#I tried this method of parameter combo for DCE on 6/2/2020, but I had issues with it

##param_combos = [] #array that will contain 1 copy of each unique parameter combo
##col_ind = np.zeros((5,1)) #column indices for folders array
###-1 is default non-filled value
##folders = -1*np.ones((5,10)) #each row is different param combo, each column is different folder for that param combo
##folders_ind = -1*np.ones((5,10)) #each row is different param combo, each column is different folder for that param combo
##
##
###loop that fills unique parameter combos array and folders (and indices) for each parameter combo
##for a in range(len(dce_ind)):
##    fold_inf = all_folders_info[dce_ind[a]]
##    paramstr = str(fold_inf.matrix) + " " + str(fold_inf.im_in_acq) + " " + str(fold_inf.pixsize0) + " " + str(fold_inf.pixsize1) + " " + str(fold_inf.orient) + " " + str(fold_inf.fatsat) + " " + str(fold_inf.tr) + " " + str(fold_inf.fa)
##
##    if (paramstr not in param_combos):
##        param_combos.append(paramstr)
##        row_ind = len(param_combos) - 1
##        folders[row_ind,int(col_ind[row_ind,0])] = dce_folders[a]
##        folders_ind[row_ind,int(col_ind[row_ind,0])] = dce_ind[a]
##        col_ind[row_ind,0] = col_ind[row_ind,0] + 1 #update column index for 2nd entry for this parameter combo
##    else:
##        for pc in range(len(param_combos)):
##            if param_combos[pc] == paramstr:
##                row_ind = pc
##                break
##            
##        folders[row_ind,int(col_ind[row_ind,0])] = dce_folders[a]
##        folders_ind[row_ind,int(col_ind[row_ind,0])] = dce_ind[a]
##        col_ind[row_ind,0] = col_ind[row_ind,0] + 1 #update column index for next entry for this parameter combo
##
###remove rows that are non-filled
##folders = folders[0:len(param_combos)-1,:]
##folders_ind = folders_ind[0:len(param_combos)-1,:]
##
##print(folders)
###iterate through these each param combo's folders to find dce_folders
##brk_out = 0
##max_folders_row = 0
##max_folders_val = 0
##
###initialize empty t1 arrays in case not filled here
##t1_folders = []
##t1_ind = []
##
###loop that decides which parameter combos folder's are DCE and which are T1
##for r in range(np.size(folders_ind,0)):
##    #return filled values of folder # and index for this parameter combo
##    curr_folders = folders[r,:]
##    curr_folders = curr_folders[curr_folders!=-1]
##    curr_folders = curr_folders.tolist()
##    
##    curr_folders_ind = folders_ind[r,:] 
##    curr_folders_ind = curr_folders_ind[curr_folders_ind!=-1]
##    curr_folders_ind = curr_folders_ind.tolist()
##
##    #keep updating which parameter combo has max # of folders
##    if len(curr_folders) > max_folders_val:
##        max_folders_row = r
##        max_folders_val = len(curr_folders)
##
##    for c in range(len(curr_folders_ind)):
##        curr_rc_ind = int(curr_folders_ind[c])
##        fold_inf = all_folders_info[curr_rc_ind]
##
##        #If PRE/POST or PRE AND POST found in series description, assume current row is for dce and all others for t1
##        if ('PRE/POST' in fold_inf.serdesc or 'PRE AND POST' in fold_inf.serdesc):
##            dce_folders = curr_folders
##            dce_ind = curr_folders_ind
##
##            folders[r,:] = -1
##            folders_vect = folders.ravel()
##            folders_vect = folders_vect[folders_vect!=-1]
##            t1_folders = folders_vect.tolist()
##            
##            folders_ind[r,:] = -1
##            folders_ind_vect = folders_ind.ravel()
##            folders_ind_vect = folders_ind_vect[folders_ind_vect!=-1]
##            t1_ind = folders_ind_vect.tolist()
##
##            brk_out = 1 #make this 1 to indicate you should break from outer loop too
##            break
##
##    if brk_out == 1:
##        break
##
###if we didn't find 'PRE/POST' or 'PRE AND POST', assume the row with most elements is dce parameter combo and everything else is t1
##if brk_out == 0:
##    dce_folders = folders[max_folders_row,:]
##    dce_folders = dce_folders[dce_folders!=0]
##    dce_folders = dce_folders.tolist()
##    
##    dce_ind = folders_ind[max_folders_row,:]
##    dce_ind = dce_ind[dce_ind!=0]
##    dce_ind = dce_ind.tolist()
##
##    folders[max_folders_row,:] = -1
##    folders_vect = folders.ravel()
##    folders_vect = folders_vect[folders_vect!=-1]
##    t1_folders = folders_vect.tolist()
##            
##    folders_ind[max_folders_row,:] = -1
##    folders_ind_vect = folders_ind.ravel()
##    folders_ind_vect = folders_ind_vect[folders_ind_vect!=-1]
##    t1_ind = folders_ind_vect.tolist()

