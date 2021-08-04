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
#Code that runs functions for DCE folder and early & late
#post-contrast identification.
#7/27/2021: Start with Plugins folder import

import Breast_DCEMRI_FTV_plugins1
from Breast_DCEMRI_FTV_plugins1 import Get_header_info_all_manufacturer
from Breast_DCEMRI_FTV_plugins1 import chooseEarly_LateByManufacturer
from Breast_DCEMRI_FTV_plugins1 import gzip_gunzip_pyfuncs
import os
import pydicom
import dicom
import numpy as np
import slicer
import qt
import sys


def runExamIdentAndTiming(exampath,dce_folders_manual,dce_ind_manual,earlyadd,lateadd):


    #First, fill header info structures for all DICOM folders in exam
    img_folders, all_folders_info = Get_header_info_all_manufacturer.fillExamFolderInfoStructures(exampath)

    #Then, run the appropriate manufacturer-specific folder identification and timing code
    folder1info = all_folders_info[0]

    #Workaround in case manufacturer header field is blank in the first
    #folder info structure
    folder = 0
    while(folder1info.manufacturer == ''):
        folder = folder + 1
        folder1info = all_folders_info[folder]

    if ('GE' in folder1info.manufacturer):
        print("---------- GE exam identification ----------")

        #Edit 2/12/2021: Only run DCE folder ident if dce_folders_manual (list of user's manual DCE folder selections) is empty
        if(len(dce_folders_manual) == 0):
            all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind = Get_header_info_all_manufacturer.ge_folder_lbl(img_folders, all_folders_info)
        else:
            dce_folders = dce_folders_manual
            dce_ind = dce_ind_manual

        #4/26/21: If you did auto DCE ID and it returned 0 DCE folders, exit this function
        if(len(dce_folders_manual) == 0 and len(dce_folders) == 0):
            slicer.util.confirmOkCancelDisplay("Automatic series ID failed. Please try manual series ID for this exam","Series ID Failed")




        #Edit 11/25/2020: commenting out 11/24 edit because I noticed that temporal resolution is on the order of ~95 seconds
        #Edit 11/24/2020: Add fix to change from 2 folders to 1 if gap between 1st and 2nd phases of multivolume folder is > 10 seconds
##        if(len(dce_folders) == 2):
##            if('gunzipped' in exampath):
##                orig_exampath = exampath[:-10] #remove \gunzipped from the end to get original exampath
##                gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[0]))
##                gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[1]))
##
##            if(folder1info.totslices/folder1info.im_in_acq > 1):
##                multiinfo = folder1info
##                multifolder = img_folders[dce_ind[0]]
##                multiind = dce_ind[0]
##            else:
##                folder2info = all_folders_info[1]
##                if(folder2info.totslices/folder2info.im_in_acq > 1):
##                    multiinfo = folder2info
##                    multifolder = img_folders[dce_ind[1]]
##                    multiind = dce_ind[1]
##
##            multipath = os.path.join(exampath,str(multifolder))
##            multifiles = os.listdir(multipath)
##            multictimes = []
##            multitrigtimes = []
##            for mf in range(len(multifiles)):
##                curr_file = multifiles[mf]
##                curr_path = os.path.join(multipath,curr_file)
##                try:
##                    curr_hdr = pydicom.dcmread(curr_path,stop_before_pixels = True)
##                except:
##                    curr_hdr = dicom.read_file(curr_path)
##                curr_ctime = curr_hdr.ContentTime
##                curr_ctimesec = 3600*int(curr_ctime[0:2]) + 60*int(curr_ctime[2:4]) + float(curr_ctime[4:])
##                multictimes.append(curr_ctimesec)
##                multitrigtimes.append(curr_hdr.TriggerTime/1000) #convert to seconds
##
##            multictimes = np.array(multictimes)
##            multictimes = np.unique(multictimes)
##
##            #For some exams, like UKCC, use trigger time in place of content time
##            if(len(multictimes) == 1):
##                multitrigtimes = np.array(multitrigtimes)
##                multictimes = np.unique(multitrigtimes)
##
##            p1p2diff = float(multictimes[1]) - float(multictimes[0])
##            if(p1p2diff > 10):
##                print("There is pause of > 10 sec between phases 1 and 2 of multivolume folder. Therefore, this is the only DCE folder.")
##                dce_folders = [multifolder]
##                dce_ind = [multiind]
        #Edit 11/24/2020 ends here


        #folder1info = all_folders_info[0]
        manufacturer = folder1info.manufacturer
        try:
            tempres, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss = chooseEarly_LateByManufacturer.chooseEarlyLateGE_Siemens(exampath,dce_folders,manufacturer,earlyadd,lateadd)
        except:
            if(len(dce_folders_manual) == 0):
                slicer.util.confirmOkCancelDisplay("Error occurred during early/late timing ID. Please try manual series ID for this exam.","Early/Late Timing ID Error")
            else:
                slicer.util.confirmOkCancelDisplay("Error occurred during early/late timing ID. If your manual series selection is correct, Slicer FTV cannot process this exam.","Early/Late Timing ID Error")

        if len(dce_folders) > 2:
            print( "Pre-contrast image is in folder #" + str(dce_folders[0]) )
            print( "Early post-contrast image is in folder #" + str(dce_folders[earlyPostContrastNum]) )
            print( "Late post-contrast image is in folder #" + str(dce_folders[latePostContrastNum]) )

        if len(dce_folders) == 1:
            print("All DCE images are in folder #" + str(dce_folders[0]))
            print( "Early post-contrast image is post-contrast #" + str(earlyPostContrastNum) )
            print( "Late post-contrast image is post-contrast #" + str(latePostContrastNum) )

        if len(dce_folders) == 2:
            print( "Pre-contrast image is in folder #" + str(dce_folders[0]) )
            print("All post-contrast images are in folder #" + str(dce_folders[1]))
            print( "Early post-contrast image is post-contrast #" + str(earlyPostContrastNum) )
            print( "Late post-contrast image is post-contrast #" + str(latePostContrastNum) )



    if ('SIEMENS' in folder1info.manufacturer or 'Siemens' in folder1info.manufacturer):
        print("---------- Siemens exam identification ----------")
        #Edit 2/12/2021: Only run DCE folder ident if dce_folders_manual (list of user's manual DCE folder selections) is empty
        if(len(dce_folders_manual) == 0):
            all_folders_info, dce_folders, dce_ind, t1_folders, t1_ind, dce_sub_folders, dce_sub_ind = Get_header_info_all_manufacturer.siemens_folder_lbl(img_folders, all_folders_info)
        else:
            dce_folders = dce_folders_manual
            dce_ind = dce_ind_manual

        #4/26/21: If you did auto DCE ID and it returned 0 DCE folders, exit this function
        if(len(dce_folders_manual) == 0 and len(dce_folders) == 0):
            slicer.util.confirmOkCancelDisplay("Automatic series ID failed. Please try manual series ID for this exam","Series ID Failed")

        #folder1info = all_folders_info[0]
        manufacturer = folder1info.manufacturer

        try:
            tempres, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss = chooseEarly_LateByManufacturer.chooseEarlyLateGE_Siemens(exampath,dce_folders,manufacturer,earlyadd,lateadd)
        except:
            if(len(dce_folders_manual) == 0):
                slicer.util.confirmOkCancelDisplay("Error occurred during early/late timing ID. Please try manual series ID for this exam.","Early/Late Timing ID Error")
            else:
                slicer.util.confirmOkCancelDisplay("Error occurred during early/late timing ID. If your manual series selection is correct, Slicer FTV cannot process this exam.","Early/Late Timing ID Error")

        if len(dce_folders) > 1:
            print( "Pre-contrast image is in folder #" + str(dce_folders[0]) )
            print( "Early post-contrast image is in folder #" + str(dce_folders[earlyPostContrastNum]) )
            print( "Late post-contrast image is in folder #" + str(dce_folders[latePostContrastNum]) )

        if len(dce_folders) == 1:
            print("All DCE images are in folder #" + str(dce_folders[0]))
            print( "Early post-contrast image is post-contrast #" + str(earlyPostContrastNum) )
            print( "Late post-contrast image is post-contrast #" + str(latePostContrastNum) )


    if ('PHILIPS' in folder1info.manufacturer or 'Philips' in folder1info.manufacturer):
        print("---------- Philips exam identification ----------")
        #Edit 2/12/2021: Only run DCE folder ident if dce_folders_manual (list of user's manual DCE folder selections) is empty
        if(len(dce_folders_manual) == 0):
            print("Using auto series ID")
            all_folders_info, dce_folders, dce_ind = Get_header_info_all_manufacturer.philips_folder_lbl(img_folders, all_folders_info)
        else:
            ("Using manual series ID")
            dce_folders = dce_folders_manual
            dce_ind = dce_ind_manual

        #4/26/21: If you did auto DCE ID and it returned 0 DCE folders, exit this function
        if(len(dce_folders_manual) == 0 and len(dce_folders) == 0):
            slicer.util.confirmOkCancelDisplay("Automatic series ID failed. Please try manual series ID for this exam","Series ID Failed")

        print("DCE folders are")
        print(dce_folders)
        print("Running timing ID")

        try:
            tempres, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss = chooseEarly_LateByManufacturer.chooseEarlyLatePhilips(exampath,dce_folders,earlyadd,lateadd)
        except:
            if(len(dce_folders_manual) == 0):
                slicer.util.confirmOkCancelDisplay("Error occurred during early/late timing ID. Please try manual series ID for this exam.","Early/Late Timing ID Error")
            else:
                slicer.util.confirmOkCancelDisplay("Error occurred during early/late timing ID. If your manual series selection is correct, Slicer FTV cannot process this exam.","Early/Late Timing ID Error")

        if(len(dce_folders) == 1):
            print("All DCE images are in folder #" + str(dce_folders[0]))
            print( "Early post-contrast image is post-contrast #" + str(earlyPostContrastNum) )
            print( "Late post-contrast image is post-contrast #" + str(latePostContrastNum) )
        if(len(dce_folders) == 2):
            print( "Pre-contrast image is in folder #" + str(dce_folders[0]) )
            print("All post-contrast images are in folder #" + str(dce_folders[1]))
            print( "Early post-contrast image is post-contrast #" + str(earlyPostContrastNum) )
            print( "Late post-contrast image is post-contrast #" + str(latePostContrastNum) )


    return tempres, all_folders_info, dce_folders, dce_ind, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss

