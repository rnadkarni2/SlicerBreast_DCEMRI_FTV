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
#Code that generates PE map, SER map,
#and tumor mask that is used to compute
#functional tumor volume (FTV).
import Breast_DCEMRI_FTV_plugins2
from Breast_DCEMRI_FTV_plugins2 import read_DCE_images_to_numpy
from Breast_DCEMRI_FTV_plugins2 import multivolume_folder_sort

import numpy as np
from scipy import signal
import nibabel as nib
import vtk
import os

def makeFTVMaps(exampath, manufacturer, dce_folders, roicenter,roiradius,omitcenters,omitradii, earlyPostContrastNum, latePostContrastNum,pct,pethresh,minconnpix):
    #Edit 10/6/2020: Make background threshold % (pct), pethresh, and minconnpix inputs to the function, rather than fixed values.
    #Doing this allows user to adjust thresholds in Slicer FTV module 2

    #Edit 6/9/2020: Even for GE or Siemens, use Philips method of loading images into numpy array if all DCE images are in same folder

    print("Running make FTV maps")
    print("dce folders")
    print(dce_folders)
    #Loading pre-contrast image into numpy array
    if ('PHILIPS' in manufacturer or 'Philips' in manufacturer or len(dce_folders) <= 2):
        #For UKCC, which has 2 DCE folders
        if(len(dce_folders) == 2):
            dcepath = os.path.join(exampath,str(dce_folders[1]))
        else:
            dcepath = os.path.join(exampath,str(dce_folders[0]))
        print("dce path is")
        print(dcepath)
        fsort, numtemp, nslice, ctime_forphases, ttime_forphases, atime_forphases = multivolume_folder_sort.sortDicoms(dcepath)

        #For UKCC, which has 2 DCE folders
        if(len(dce_folders) == 2):
            apath = os.path.join(exampath,str(dce_folders[0]))
            print("pre-contrast path is")
            print(apath)
            m,a = read_DCE_images_to_numpy.readInputToNumpy(apath)
        else:
            m,a = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,0)

    else:
        apath = os.path.join(exampath,str(dce_folders[0])) #5/8/2020: Edited to go with identify_dce_folders
        m,a = read_DCE_images_to_numpy.readInputToNumpy(apath)


    #Loading early post-contrast image into numpy array
    print("-----EARLY POST-CONTRAST IMAGE-----")
    if ('PHILIPS' in manufacturer or 'Philips' in manufacturer or len(dce_folders) <= 2):
      m,b = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,earlyPostContrastNum)
    else:
      m,b = read_DCE_images_to_numpy.earlyOrLateImgSelect(earlyPostContrastNum,dce_folders,exampath) #5/8/2020: Edited to go with identify_dce_folders

    #function used to add an ROI or omit region to a voi_mask numpy array
    def addToVOIMask(voi_mask,center,radius,maskval):
      #maskval: 1 if ROI, 0 if omit region

      #tranpose center and radius to avoid indexing errors
      center = np.transpose(center)
      radius = np.transpose(radius)

      #first, get x,y,z range of region from its center and radius
      xs = int( round( float(center[0]) - float(radius[0]) ) )
      xf = int( round( float(center[0]) + float(radius[0]) ) )

      ys = int( round( float(center[1]) - float(radius[1]) ) )
      yf = int(round(float(center[1]) + float(radius[1])))

      zs = int( round( float(center[2]) - float(radius[2]) ) )
      zf = int( round( float(center[2]) + float(radius[2]) ) )

      #because it's a numpy array, it must have dimension order z,y,x in order to match orientations with input image
      #4/27/21: Add the +1 to each dimension to compensate for
      #Python array indexing rules
      voi_mask[xs:xf+1,ys:yf+1,zs:zf+1] = maskval
      return xs,xf,ys,yf,zs,zf,voi_mask

    #same dimension order as other images -- x,y,z
    voi_mask = np.zeros(a.shape)

    #add roi to voi mask
    xs,xf,ys,yf,zs,zf,voi_mask = addToVOIMask(voi_mask,roicenter,roiradius,1)

    #add omit regions to voi mask
    if(omitcenters[0,0] != -1):
        for i in range(omitcenters.shape[0]):
            current_omit_center = omitcenters[i,:]
            current_omit_radius = omitradii[i,:]

            #for now, structured to only keep coords of last omit
            oxs,oxf,oys,oyf,ozs,ozf,voi_mask = addToVOIMask(voi_mask,current_omit_center,current_omit_radius,0)

    #convert precontrast image's numpy array into VOI cropped, masked version
    #Edit 5/20/2020: correct dimension order is x,y,z--verified by correct cropped outputs
    #with this dimension order in DCE_TumorMapProcess
    amasked = a*voi_mask
    amasked = amasked[xs:xf+1,ys:yf+1,zs:zf+1]

    #4/27/21: Add the +1 to each dimension to compensate for
    #Python array indexing rules
    
    #use roi cropped omit region zeroed version of pre-contrast image to compute pre-contrast threshold
    #pct = 0.6 #minimum percent of max (or in this case 95%ile) value in VOI used to define pre-contrast threshold
    pre_thresh = pct*np.percentile(amasked,95) #0.6 x 95th percentile value in ROI
    #don't use max instead of 95%ile -- max gives extremely high value of pre_thresh and makes tumor volume 0
    print("BKG threshold value")
    print(pre_thresh)

    pe = 100*(b-a)/a
    pe = np.where(a>=pre_thresh,pe,0) #set PE to 0 where pre-contrast is less than BKG thresh

    #Loading late post-contrast image into numpy array
    print("-----LATE POST-CONTRAST IMAGE-----")
    if ('PHILIPS' in manufacturer or 'Philips' in manufacturer or len(dce_folders) <= 2):
      m,c = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,latePostContrastNum)
    else:
      m,c = read_DCE_images_to_numpy.earlyOrLateImgSelect(latePostContrastNum,dce_folders,exampath) #5/8/2020: Edited to go with identify_dce_folders

    ser = (b-a)/(c-a)
    ser = np.where(a>=pre_thresh,ser,0)

    #Creating tumor mask
    br_mask = (a>=pre_thresh)

    pe_mask = (pe>=pethresh)

    #5/8/21: Using this as convolution input because doing so reduced
    #Slicer FTV bias compared to brtool
    br_pe_mask = br_mask*pe_mask

    kernel = np.ones((3,3,3))
    kernel[1,1,1] = 100

    convbrmask = signal.convolve(br_pe_mask,kernel,mode='same')

    #add affine matrix IJKToRAS to numpy array so that it can be added to RGB image
    aff_mat_RAS = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            aff_mat_RAS[i,j] = m.GetElement(i,j)

    connpix_mask = (convbrmask>=(100+minconnpix))

    ser_mask = (ser >= 0) #4/28/21: Explicitly use SER>=0 for FTV definition

    tumor_mask = br_pe_mask*connpix_mask*ser_mask #tumor segment #Edit 4/28/21: multiply by ser_mask too
    tumor_mask = tumor_mask.astype('float64') #Convert from bool to numeric to prevent NIFTI error

    print("Done running make FTV maps")

    return a,b,c,pe,ser,tumor_mask,voi_mask,zs,zf,ys,yf,xs,xf,pct,pre_thresh,pethresh,minconnpix


