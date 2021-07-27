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

def makeFTVMaps(exampath, manufacturer, dce_folders, roicenter,roiradius,omitcenters,omitradii, earlyPostContrastNum, latePostContrastNum,pct,pethresh,minconnpix):
    #Edit 10/6/2020: Make background threshold % (pct), pethresh, and minconnpix inputs to the function, rather than fixed values.
    #Doing this allows user to adjust thresholds in Slicer FTV module 2
    

    #Call function that tells which images are early and late post-contrast
    #Edit 5/8/2020: Commenting out line below because need to make timing code that uses outputs from identify_dce_folders functions
    #studydate, earlyPostContrastNum,latePostContrastNum,prefoldernum, earlydiffmm, earlydiffss, latediffmm, latediffss = choose_early_late_imgs_func.chooseEarlyLate(exampath)
    #Edit 6/9/2020: Even for GE or Siemens, use Philips method of loading images into numpy array if all DCE images are in same folder

    print("Running make FTV maps")
    print("dce folders")
    print(dce_folders)
    #Loading pre-contrast image into numpy array
    if ('PHILIPS' in manufacturer or 'Philips' in manufacturer or len(dce_folders) <= 2):
        #For UKCC, which has 2 DCE folders
        if(len(dce_folders) == 2):
            dcepath = exampath + "\\" + str(dce_folders[1])
        else:
            dcepath = exampath + "\\" + str(dce_folders[0])
        print("dce path is")
        print(dcepath)
        fsort, numtemp, nslice, ctime_forphases, ttime_forphases, atime_forphases = multivolume_folder_sort.sortDicoms(dcepath)

        #For UKCC, which has 2 DCE folders
        if(len(dce_folders) == 2):
            apath = exampath + "\\" + str(dce_folders[0])
            print("pre-contrast path is")
            print(apath)
            m,a = read_DCE_images_to_numpy.readInputToNumpy(apath)
        else:
            m,a = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,0)
            
    else:
        apath = exampath + "\\" + str(dce_folders[0]) #5/8/2020: Edited to go with identify_dce_folders
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
    #voi_mask = voi_mask.astype('float64')

    #add roi to voi mask
    xs,xf,ys,yf,zs,zf,voi_mask = addToVOIMask(voi_mask,roicenter,roiradius,1)

    #add omit regions to voi mask
    if(omitcenters[0,0] != -1):
        for i in range(omitcenters.shape[0]):
            current_omit_center = omitcenters[i,:]
            current_omit_radius = omitradii[i,:]

            #for now, structured to only keep coords of last omit
            oxs,oxf,oys,oyf,ozs,ozf,voi_mask = addToVOIMask(voi_mask,current_omit_center,current_omit_radius,0)


##    print("# of voxels in VOI mask")
##    print(np.sum(voi_mask))

    #convert precontrast image's numpy array into VOI cropped, masked version
    #Edit 5/20/2020: correct dimension order is x,y,z--verified by correct cropped outputs
    #with this dimension order in DCE_TumorMapProcess
##    try:
    amasked = a*voi_mask
    amasked = amasked[xs:xf+1,ys:yf+1,zs:zf+1]
##    except:
##        amasked = np.zeros((xf+1-xs,yf+1-ys,zf+1-zs))
##        a2 = a.astype('int8')
##        a2 = a2[xs:xf+1,ys:yf+1,zs:zf+1]
##        voi_mask2 = voi_mask.astype('int8')
##        voi_mask2 = voi_mask2[xs:xf+1,ys:yf+1,zs:zf+1]
##        amasked = a2*voi_mask
##        amasked = amasked.astype('float64')
        
    #4/27/21: Add the +1 to each dimension to compensate for
    #Python array indexing rules
    #amasked = amasked.astype('int16') #4/29/21: seeing if data type affects results

##    print("amasked shape")
##    print(amasked.shape)

    #use roi cropped omit region zeroed version of pre-contrast image to compute pre-contrast threshold
    #pct = 0.6 #minimum percent of max (or in this case 95%ile) value in VOI used to define pre-contrast threshold
    
    pre_thresh = pct*np.percentile(amasked,95) #0.6 x 95th percentile value in ROI
    #don't use max instead of 95%ile -- max gives extremely high value of pre_thresh and makes tumor volume 0
    print("BKG threshold value")
    print(pre_thresh)


    #7/13/2021: Just crop the images to include only VOI when memory error occurs.
    #Trying other data types doesn't seem to work that well.
##    try:
    pe = 100*(b-a)/a
##    except:
##        acrop = a[xs:xf+1,ys:yf+1,zs:zf+1]
##        bcrop = b[xs:xf+1,ys:yf+1,zs:zf+1]
##        pe = 100*(bcrop-acrop)/(acrop)

    #Undid the following -> 4/30/21: See if commenting out this step increases similarity to brtool copy
    pe = np.where(a>=pre_thresh,pe,0) #set PE to 0 where pre-contrast is less than BKG thresh
    #pe = pe.astype('int16') #4/29/21: seeing if data type affects results
    #pe = pe.astype('int16') #Edit 5/12/2020: only need 64 bit for calculation. convert back to 16 bit to decrease storage space


    #Loading late post-contrast image into numpy array
    print("-----LATE POST-CONTRAST IMAGE-----")
    if ('PHILIPS' in manufacturer or 'Philips' in manufacturer or len(dce_folders) <= 2):
      m,c = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,latePostContrastNum)
    else:
      m,c = read_DCE_images_to_numpy.earlyOrLateImgSelect(latePostContrastNum,dce_folders,exampath) #5/8/2020: Edited to go with identify_dce_folders

    #7/13/2021: Just crop the images to include only VOI when memory error occurs.
    #Trying other data types doesn't seem to work that well.
##    try:
    ser = (b-a)/(c-a)
    ser = np.where(a>=pre_thresh,ser,0)
##    except:
##        acrop = a[xs:xf+1,ys:yf+1,zs:zf+1]
##        bcrop = b[xs:xf+1,ys:yf+1,zs:zf+1]
##        ccrop = c[xs:xf+1,ys:yf+1,zs:zf+1]
##        ser = (bcrop-acrop)/(ccrop-acrop)
##        ser = np.where(acrop>=pre_thresh,ser,0)

        
    #4/30/21: ser has to be float64 (cannot be int) to give normal results for
    #SER color breakdown
    #ser = ser.astype('float64') #4/29/21: seeing if data type affects results

    #make tumor mask
    #percent enhancement threshold
    #pethresh = 70
    #minimum connected pixels threshold
    #minconnpix = 3
    
    #Creating tumor mask
##    try:
    br_mask = (a>=pre_thresh)
##    except:
##        a2 = a.astype('int8')
##        br_mask = (a2>=pre_thresh)
    #br_mask = br_mask.astype('int16') #4/29/21: seeing if data type affects results
##    print("# of voxels in VOI that are above BKG threshold")
##    print(np.sum(br_mask*voi_mask))

    #7/13/2021: Need to make sure that if one of the maps is cropped,
    #all are cropped.
##    if(br_mask.shape[0] > ser.shape[0] or br_mask.shape[0] > pe.shape[0]):
##        #first, crop the VOI mask to only include VOI
##        #(now it's all 1's unless there's an omit region).
##        voi_mask = voi_mask[xs:xf+1,ys:yf+1,zs:zf+1]
##        #then, crop the br_mask to only include VOI
##        br_mask = br_mask[xs:xf+1,ys:yf+1,zs:zf+1]
##
##        #if SER is cropped but PE isn't, crop PE
##        if(pe.shape[0] > ser.shape[0]):
##            pe = pe[xs:xf+1,ys:yf+1,zs:zf+1]
##            
##        #if PE is cropped but SER isn't, crop SER
##        if(ser.shape[0] > pe.shape[0]):
##            ser = ser[xs:xf+1,ys:yf+1,zs:zf+1]
    
    pe_mask = (pe>=pethresh)

    #5/8/21: Using this as convolution input because doing so reduced
    #Slicer FTV bias compared to brtool
    br_pe_mask = br_mask*pe_mask

    kernel = np.ones((3,3,3))
    kernel[1,1,1] = 100

##    try:
    convbrmask = signal.convolve(br_pe_mask,kernel,mode='same')
##    except:
##        br_pe_mask = br_pe_mask[xs:xf+1,ys:yf+1,zs:zf+1]
##        convbrmask = signal.convolve(br_pe_mask,kernel,mode='same')
##        voi_mask = voi_mask[xs:xf+1,ys:yf+1,zs:zf+1] 
##        ser = ser[xs:xf+1,ys:yf+1,zs:zf+1]

    #Temporary code for saving convbrmask to nii file

    #add affine matrix IJKToRAS to numpy array so that it can be added to RGB image
    aff_mat_RAS = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            aff_mat_RAS[i,j] = m.GetElement(i,j)

    #save convbrmask to nifti with affine matrix extracted from the loop above
    #convbrmask_img = nib.Nifti1Image(convbrmask, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
    #convbrmask_img_savename = r"C:\Users\rnadkarni\convolutions\conv_Slicer.nii"
    #nib.save(convbrmask_img,convbrmask_img_savename)
    
    connpix_mask = (convbrmask>=(100+minconnpix))
    #connpix_mask = connpix_mask.astype('int16') #4/29/21: seeing if data type affects results
    #connpix_mask= connpix_mask.astype('float64') #because bool can't be saved to .nii file
##    print("# of voxels in VOI that have at least 4 connected neighbors in BKG mask")
##    print(np.sum(connpix_mask*voi_mask))
##
##    print("# of voxels in VOI that are above PE threshold and have at least 4 neighbors")
##    print(np.sum(pe_mask*connpix_mask*voi_mask))

    ser_mask = (ser >= 0) #4/28/21: Explicitly use SER>=0 for FTV definition

    #7/13/2021: Need to do this differently if you had to crop
    #to save memory earlier.
##    if(a.shape[0] > voi_mask.shape[0]):
##        tumor_mask = np.zeros(a.shape)
##        tumor_mask[xs:xf+1,ys:yf+1,zs:zf+1] = br_pe_mask*connpix_mask*ser_mask
##    else:
    tumor_mask = br_pe_mask*connpix_mask*ser_mask #tumor segment #Edit 4/28/21: multiply by ser_mask too
    tumor_mask = tumor_mask.astype('float64') #Convert from bool to numeric to prevent NIFTI error
        
    print("Done running make FTV maps")


    #4/27/21: save key maps to .nii files to help David debug brtool - Slicer discrepancy
    #for 19815 v20
##    voi_mask_img = nib.Nifti1Image(voi_mask, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    voi_mask_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\VOImask_Slicer.nii"
##    nib.save(voi_mask_img,voi_mask_img_savename)
##    
##    pre_img = nib.Nifti1Image(a, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    pre_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\precontrast_Slicer.nii"
##    nib.save(pre_img,pre_img_savename)
##
##    pecrop = pe[0:256,127:383,2:162] #for 19768 v20 only
##    pe_img = nib.Nifti1Image(pecrop, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    pe_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\PE_Slicer.nii"
##    nib.save(pe_img,pe_img_savename)
##
##    #4/30/21: Save pemask too to debug brtool discrepancies.
##    pe_mask = pe_mask.astype('int16') #because bool can't be saved to .nii
##    pemask_crop = pe_mask[0:256,127:383,2:162] #for 19768 v20 only
##    pemask_img = nib.Nifti1Image(pemask_crop, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    pemask_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\PEmask_Slicer.nii"
##    nib.save(pemask_img,pemask_img_savename)
##
##
##    ser_img = nib.Nifti1Image(ser, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    ser_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\SER_Slicer.nii"
##    nib.save(ser_img,ser_img_savename)
##
##    br_mask_crop = br_mask[0:256,127:383,2:162] #for 19768 v20 only
##    br_mask_img = nib.Nifti1Image(br_mask_crop, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    br_mask_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\bkg_mask_Slicer.nii"
##    nib.save(br_mask_img,br_mask_img_savename)
##
##    connpix_mask_crop = connpix_mask[0:256,127:383,2:162] #for 19768 v20 only
##    connpix_mask_img = nib.Nifti1Image(connpix_mask_crop, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    connpix_mask_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\mnc_mask_Slicer.nii"
##    nib.save(connpix_mask_img,connpix_mask_img_savename)
##
##    tumor_mask_img = nib.Nifti1Image(tumor_mask, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
##    tumor_mask_img_savename = r"C:\Users\rnadkarni\SlcExt\FTV_process_complete\DCE_TumorMapProcess\FTV_processing_images\tumormask_Slicer.nii"
##    nib.save(tumor_mask_img,tumor_mask_img_savename)


    return a,b,c,pe,ser,tumor_mask,voi_mask,zs,zf,ys,yf,xs,xf,pct,pre_thresh,pethresh,minconnpix


