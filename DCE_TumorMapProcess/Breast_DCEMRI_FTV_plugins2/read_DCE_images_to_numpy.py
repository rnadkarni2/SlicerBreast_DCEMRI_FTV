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
#Script that contains all functions for reading
#DICOM series from a path.

import numpy as np
import pydicom
import dicom
import re
import os
import slicer
import vtk

#Philips correct z-order of slices
#Edit 6/10/2020: By checking both OHSC and UChic exams from Philips, I realized that the routine from
#multivolume_folder_sort always leads to the correct orientation in the report, regardless of if the
#filenames are in reverse order or not -- This is for Philips only


#GE & Siemens correct z-order of slices
#Update--It appears that regular, increasing order of Slice Location works for RMH exams, but not for UCSD exams
#Maybe just use the slice names for sorting again, but only for GE & Siemens?
    #-Update: Order of increasing slice name appears to work well for GE and Siemens exams --
    #gives correct orientation and FTV in report and ROI in Slicer matches well with z-range of tumor

#4/3/2020: New version of readInputToNumpy using SimpleITK
#Edit 6/16/2020: Incorporate Andrey's suggested method of reading DICOM series into Slicer
def readInputToNumpy(path):
    #Edit 6/16/2020: os.listdir and adding full path may be faster
    dicom_names = os.listdir(path)
    dicom_names = addFullPathToFileList(path,dicom_names) #call function for converting list of DCM filenames to list of DCM file names with path included
        
    plugin = slicer.modules.dicomPlugins['DICOMScalarVolumePlugin']()
    loadables = plugin.examine([dicom_names])
    volume = plugin.load(loadables[0])

    #print(loadables)
    #print(volume)

    #store RASToIJKMatrix in m so that you can use this to correct output images' orientations
    m = vtk.vtkMatrix4x4()
    volume.GetRASToIJKMatrix(m)

    npimg = slicer.util.arrayFromVolume(volume)
    print("image dimensions")
    print(np.shape(npimg))
    print("image min and max after arrayFromVolume")
    print(np.amin(npimg))
    print(np.amax(npimg))
    slicer.mrmlScene.RemoveNode(volume) #this is how to prevent "unnamed volume ..." from staying in slicer window
    npimg = np.transpose(npimg,(2,1,0)) #edit 4/27/2020: using dimension order x,y,z gives same orientation as input DICOM (original numpy one is z,y,x)
    #npimg = np.flip(npimg,2) #apparently slicer method of reading DICOM series always reverses slice order from what we want, so add this flip
    npimg = npimg.astype('float64') #convert to float64 to avoid dynamic range issues when calculating PE and SER
    print("image min and max after float64 conversion")
    print(np.amin(npimg))
    print(np.amax(npimg))
    return m, npimg


#Old sitk code for reading file
##    dicom_names = sortByName(dicom_names)
##    reader.SetFileNames(dicom_names)
##    image = reader.Execute()
##    npimg = sitk.GetArrayFromImage(image)
##    npimg = np.transpose(npimg,(2,1,0)) #edit 4/27/2020: using dimension order x,y,z gives same orientation as input DICOM (original numpy one is z,y,x)
##    #Make sure that when you take ROI from numpy array, you are following the x,y,z dimension order
##    npimg = npimg.astype('float64') #convert to float64 to avoid dynamic range issues when calculating PE and SER


#5/8/2020: Edited to go with identify_dce_folders
def earlyOrLateImgSelect(postContrastNum,dce_folders,exampath):
    #Function called to select early & late post-contrast images based on timing
    #postContrastNum is number of post-contrast image eg if folder is 1002, postContrastNum = 2
    #7/6/2021: Add exception cases for DCE folders that have other
    #characters besides numbers, like Duke TCIA ones.
    try:
        foldernum = int(dce_folders[postContrastNum])
    except:
        foldernum = str(dce_folders[postContrastNum])
        
    print("Reading post-contrast image # " + str(postContrastNum) + " from folder " + str(foldernum))
    path = os.path.join(exampath,str(foldernum))
    a = readInputToNumpy(path)
    return a

#For this function, need loop to concatenate full path to filename at each iteration
def readPhilipsImageToNumpy(exampath,dce_folders,fsort,postContrastNum):

    #If 2 dce folders, need to subtract 1 to obtain correct row of fsort because it doesn't have a row for precontrast
    if len(dce_folders) == 2:
        postContrastNum = postContrastNum - 1
        dcepath = os.path.join(exampath,str(dce_folders[1]))
        print("post-contrast dce path is")
        print(dcepath)
    else:
        dcepath = os.path.join(exampath,str(dce_folders[0]))
    
    dicom_names = fsort[postContrastNum][:]
    dicom_names = addFullPathToFileList(dcepath,dicom_names) #call function for converting list of DCM filenames to list of DCM file names with path included

    #7/20/2021: If image dimensions > 700 x 700 x 200, cut off
    #20 slices from each side. Hopefully this prevents numpy
    #memory errors when computing FTV.
    try:
        hdr0 = pydicom.dcmread(str(dicom_names[0]),stop_before_pixels = True)
        nrow = hdr0[0x28,0x10].value
        if(nrow > 700 and len(dicom_names) > 200):
            dicom_names = dicom_names[30:-30]
    except:
        print("Not checking image size, therefore will not remove slices.")

    plugin = slicer.modules.dicomPlugins['DICOMScalarVolumePlugin']()
    loadables = plugin.examine([dicom_names])
    volume = plugin.load(loadables[0])

    #print(loadables)
    #print(volume)

    #store RASToIJKMatrix in m so that you can use this to correct output images' orientations
    m = vtk.vtkMatrix4x4()
    volume.GetRASToIJKMatrix(m)
    
    npimg = slicer.util.arrayFromVolume(volume)
    print("image min and max after arrayFromVolume")
    print(dcepath)
    print(np.amin(npimg))
    print(np.amax(npimg))
    slicer.mrmlScene.RemoveNode(volume) #this is how to prevent "unnamed volume ..." from staying in slicer window

    #must do this with sitk method, otherwise there are orientation and intensity differences
    #reader = sitk.ImageSeriesReader()
    #reader.SetFileNames(dicom_names)
    #image = reader.Execute()
    #npimg = sitk.GetArrayFromImage(image)
##    for i in range(len(dicom_names)):
##        curr_path = exampath + "\\" + str(int(dce_folders[0])) + "\\" + dicom_names[i]
##        try:
##          curr_img = pydicom.dcmread(curr_path)
##        except:
##          curr_img = dicom.read_file(curr_path)
##        curr_img = curr_img.pixel_array
##        if i == 0:
##          npimg = np.zeros(( len(dicom_names), curr_img.shape[0], curr_img.shape[1]))
##        npimg[i,:,:] = curr_img

    npimg = np.transpose(npimg,(2,1,0)) #edit 4/27/2020: using dimension order x,y,z gives same orientation as input DICOM (original numpy one is z,y,x)
    #npimg = np.flip(npimg,2) #apparently slicer method of reading DICOM series always reverses slice order from what we want, so add this flip
    try:
        npimg = npimg.astype('float64') #convert to float64 to avoid dynamic range issues when calculating PE and SER
    except:
        try:
            npimg = npimg.astype('float32')
        except:
            npimg = npimg.astype('int8')
            
    return m, npimg

#Wrote this function for quickly converting list of DCM filenames to list of DCM file names with path included
def addFullPathToFileList(path,files):
    for i in range(len(files)):
        filei = files[i]
        fullfilei = os.path.join(path,filei)
        files[i] = fullfilei
    return files

#Edit 6/10/2020: Introduced this function to make sure that GE and Siemens DICOM slices are always in correct order
##def sortBySliceLocation(dicom_names):
##    slclocs = []
##    for i in range(len(dicom_names)):
##        curr_path = dicom_names[i]
##        try:
##            curr_hdr = pydicom.dcmread(curr_path,stop_before_pixels = True)
##        except:
##            curr_hdr = dicom.read_file(curr_path)
##        slcloc = float(curr_hdr[0x20,0x1041].value)
##        slclocs.append(slcloc)
##    slclocs = np.array(slclocs)
##    
##    #inds are the indices to sort dicom_names in order of slice location
##    inds = slclocs.argsort()
##
##    #loop to resort dicom_names according to inds
##    dicom_names_old = dicom_names
##    dicom_names = []
##    for j in range(len(inds)):
##        dicom_names.append(dicom_names_old[int(inds[j])])
##    return dicom_names

#Edit 6/12/2020: Add exception so it can handle the kind of dicom numbering found in UCSF ISPY ID 16078
##def sortByName(dicom_names):
##    slc1 = dicom_names[0]
##    
##    slcend = dicom_names[len(dicom_names)-1]
##
##    try:
##        slc1num = int(slc1[-7:-4])
##        slcendnum = int(slcend[-7:-4])
##    except:
##        str1 = slc1[-7:-4]
##        rr1 = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str1)
##        slc1num = int(rr1[-1])
##
##        strend = slcend[-7:-4]
##        rrend = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", strend)
##        slcendnum = int(rrend[-1])
##
##                      
##    if (slc1num > slcendnum):
##        dicom_names = list(dicom_names)
##        dicom_names.reverse()
##    return dicom_names
    
