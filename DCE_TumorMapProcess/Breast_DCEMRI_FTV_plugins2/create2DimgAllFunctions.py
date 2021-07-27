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

#File that is called to create MIP image from 3d image,
#draw ROI rectangle on 2d image,
#or find slice of 3d image with maximum tumor volume

import cv2
import numpy as np

#Function that makes MIP image
def makeMIP(img3d):
    mip = np.zeros((img3d.shape[0],img3d.shape[1]))
    for y in range(img3d.shape[0]):
        for x in range(img3d.shape[1]):
            currxyallslice = img3d[y,x,:]
            maxcurrxy = np.amax(currxyallslice)
            mip[y,x] = maxcurrxy
    return mip

#Function that draws box for ROI on 2D image
#Edit 7/24/2020: Turn image into RGB image and draw orange box
#view is a char array that should have value of either 'ax' or 'sag' to specify if input image is of axial or sagittal orientation
def createImgWithROIRect(img2d,col_s,col_f,row_s,row_f,omitCount, omitradii, omitcenters,view,slc_maxA):
    img2dcp = img2d.copy() #copy so that rect doesn't stay on original image
    img2dcp = np.ascontiguousarray(img2dcp) #prevents data type errors in cv2.rectangle line
    norm_fact = 255 #float(np.percentile(img2dcp,90))
    img2dcp = img2dcp/norm_fact #normalize by 90th percentile value -- Edit: using 255 as normalization factor instead
    img2dcp_rgb = np.dstack((img2dcp,img2dcp,img2dcp)) #convert to RGB image
    img2dcp_rgb = cv2.rectangle(img2dcp_rgb,(col_s,row_s),(col_f,row_f),[255,128,0],1) #For ROI draw orange rectangle with thickness 1 on grayscale image
    #Drawing filled orange rectangles for the omit regions
    for i in range(omitCount):
        curroc = omitcenters[i,:]
        curror = omitradii[i,:]
        currstart = curroc-curror
        currend = curroc+curror
        
        #fill in orange rectangle for ith omit region (thickness -1 = filled)
        if(view == 'ax'): #if axial, use x and y values of omit
            #only add omit to image if it would actually be on the axial slice you're displaying
            if(int(currstart[2])<=slc_maxA and slc_maxA<=int(currend[2])):
                img2dcp_rgb = cv2.rectangle(img2dcp_rgb,(int(currstart[0]),int(currstart[1])),(int(currend[0]),int(currend[1])),[255,128,0],-1)
        if(view == 'sag'): #if sagittal, use y and z values of omit
            #only add omit to image if it would actually be on the sagittal slice you're displaying
            if(int(currstart[0])<=slc_maxA and slc_maxA<=int(currend[0])):
                img2dcp_rgb = cv2.rectangle(img2dcp_rgb,(int(currstart[1]),int(currstart[2])),(int(currend[1]),int(currend[2])),[255,128,0],-1)
        
    return img2dcp_rgb

#function that finds out which slice shows the greatest tumor area so that you can choose to display this slice
def chooseMaxTumorSlice(slc_s,slc_f,tumor_mask):
    slc_sums = np.zeros(((1+slc_f-slc_s),1)) #array that contains tumor area of each slice in roi slice range
    slc_maxA = 0 #slice index of max tumor area
    maxsum = 0 #max tumor area in that slice
    for slc in range(slc_s,(slc_f+1)):
        tumor_mask_slc = tumor_mask[:,:,slc]
        #update max slice index when current tumor area sum exceeds previous max value
        if(np.sum(tumor_mask_slc) > maxsum):
            maxsum = np.sum(tumor_mask_slc)
            slc_maxA = slc
    return slc_maxA
