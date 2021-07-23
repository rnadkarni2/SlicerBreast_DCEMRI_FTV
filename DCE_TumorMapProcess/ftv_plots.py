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
#Code for creating Slicer FTV report that includes FTV, ROI boundaries,
#images from the MR exam, and other relevant information such as
#imaging site, visit number, and exam date.

def createPDFreport(gzipped,path,savenamepdf,tempres,fsort,manufacturer,dce_folders,nslice,earlyPostContrastNum,latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss, preimg3d,img3d,ser,tumor_mask,voi_mask,xs,xf,ys,yf,zs,zf,omitCount,omitradii,omitcenters,pct,pre_thresh,pethresh,minconnpix,aff_mat,ijkToRASmat,nodevisstr,window,level,idstr):
    #note: although variable is called ser_colormap, may decide to use regular SER image instead so that colorbar shows true SER values



    #All of the figures (except ser colormap) are using 2nd post-contrast minus pre-contrast
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import cv2
    import numpy as np
    from skimage import color
    import pydicom
    import dicom
    import os
    import create2DimgAllFunctions
    #import parse_xml



    #Edit 8/12/2020: Use window and level to set display min and max
    dispmin = int(level)
    dispmax = int(level+window/2)
    minmaxstr = "Displaying images with min " + str(dispmin) + " and max" + str(dispmax)
    print(minmaxstr)


    #Edit 4/27/2020: Since images now use dimension order x,y,z
    #but we still want axial images to have dimension order y,x,z for this report,
    #add tranpose step here to make dimension order y,x,z
    preimg3d = np.transpose(preimg3d,(1,0,2))
    print("precontrast min and max")
    print(np.amin(preimg3d))
    print(np.amax(preimg3d))

    img3d = np.transpose(img3d,(1,0,2))
    print("early postcontrast min and max")
    print(np.amin(img3d))
    print(np.amax(img3d))


    
    ser = np.transpose(ser,(1,0,2))
    tumor_mask = np.transpose(tumor_mask,(1,0,2))
    voi_mask = np.transpose(voi_mask,(1,0,2))

    #Edit 5/11/2020: Multiply tumor mask by voi mask to prevent SER colorization in omit regions
    #7/13/2021: Only multiply by voi_mask if maps are not cropped.
##    if(voi_mask.shape[0] == tumor_mask.shape[0]):
    tumor_mask = tumor_mask*voi_mask


    
    #Method for sagittal image to make width/height match aspect ratio
    #Method should be biased towards cutting off the chest wall behind the breast rather than cutting off nipple
##    def cropSagittal(sagimg,ys,yf,aspect):
##        space_diff = aspect*sagimg.shape[0]-(yf-ys) #difference between aspect x # of rows in image and # of columns occupied by ROI in sagittal slice
##        if space_diff>0:
##            #crop img to aspect x #rows by adding 2/3 of the space difference to end of yf and 1/3 before ys
##            col_s = ys-int(space_diff/3)
##            if col_s < 0:
##                col_s = 0
##            col_f = yf+int(2*space_diff/3)
##            if col_f > (sagimg.shape[1]-1):
##                col_f = sagimg.shape[1]-1
##            sagimg = sagimg[:,col_s:col_f]
##            #incorporate l/r flip that gives desired orientation
##        sagimg = np.fliplr(sagimg)
##        return sagimg

    #Function to colorize tumor in ROI based on SER values
    #Edit 7/24/2020: Change function to reflect fact that input is now an RGB image
    def serColorize(img_wroi,tumor_mask_slc,ser_slc,row_s,row_f,col_s,col_f):
        #edit 5/19/2020: first, re-format input image by setting negative values to zero and normalizing image to max value 255
        #img_wroi = np.where(img_wroi>0,img_wroi,0) #set to 0 in voxels where image is not greater than 0

        #img_wroi = img_wroi/np.percentile(img_wroi,90) #normalize by 90th percentile value
        
        #initialize output array by making RGB version of input slice w ROI box
        #img_wroi_sercolor = np.dstack((img_wroi,img_wroi,img_wroi))
        #img_wroi_sercolor.astype('uint8')

        #edit 5/19/2020: getting rid of old normalization method
        #img_wroi_sercolor = img_wroi_sercolor/(np.percentile(img_wroi_sercolor,pct))

        img_wroi_sercolor = np.copy(img_wroi) #need to use np.copy so that input image is not SER colorized
        
        #Only doing colorization within ROI box
        for r in range(row_s+1,row_f):
            for c in range(col_s+1,col_f):
                #Only add color to voxels that have been labeled as tumor
                if tumor_mask_slc[r,c] == 1:
                    #If SER in range 0 to 0.9, make the voxel blue
                    if ser_slc[r,c] > 0 and ser_slc[r,c] <= 0.9:
                        img_wroi_sercolor[r,c] = [0,0,255]

                    #If SER in range 0.9 to 1, make the voxel purple
                    if ser_slc[r,c] > 0.9 and ser_slc[r,c] <= 1:
                        img_wroi_sercolor[r,c] = [128,0,128]

                    #If SER in range 1 to 1.3, make the voxel green
                    if ser_slc[r,c] > 1 and ser_slc[r,c] <= 1.3:
                        img_wroi_sercolor[r,c] = [0,255,0]

                    #If SER in range 1.3 to 1.75, make the voxel red
                    if ser_slc[r,c] > 1.3 and ser_slc[r,c] <= 1.75:
                        img_wroi_sercolor[r,c] = [255,0,0]

                    #If SER above 1.75, make the voxel white
                    #Edit 5/1/2020: Make last one yellow instead of white to avoid
                    #confusion with bright pixels in grayscale image and cap at max
                    #value 3 as seen in Aegis
                    if ser_slc[r,c] > 1.75 and ser_slc[r,c]<=3:
                        img_wroi_sercolor[r,c] = [255,255,0]
                        
        return img_wroi_sercolor

        

    #Function that creates a square cropped region around ROI rectangle
    #Use this to crop axial slice with SER color values
    def sqCropROIimg(img_wroi,row_s,row_f,col_s,col_f,buf,aspect,view):
        #Edit 6/15/2020: new variable view: 'axial' or 'sagittal'
        
        #buf: buffer around larger dimension of ROI

        #row_f, row_s, col_f, col_s are the ROI rectangle bounds

        #If the ROI has aspect*(#rows)>(#columns), the goal of this
        #function is to to make the cropped
        #version of the whole image have height of (row_f-row_s) + 2*buf
        #and width of aspect*((row_f-row_s) + 2*buf)

        #If the ROI has aspect*(#rows)<(#columns), the goal is to make the
        #cropped version of the whole image have width of (col_f-col_s)+2*buf
        #and height of ((col_f-col_s)+2*buf)/aspect

        if aspect*(row_f-row_s) > (col_f-col_s):
            croprow_s = row_s-buf
            croprow_f = row_f+buf

            spacediff = aspect*(2*buf+(croprow_f-croprow_s))-(col_f-col_s)

            cropcol_s = col_s - int(0.5*spacediff)
            cropcol_f = col_f + int(0.5*spacediff)
        else:
            cropcol_s = col_s-buf
            cropcol_f = col_f+buf

            #Edit 6/15/2020: if sagittal, only crop 2*buf around row_s and row_f. For axial, keep cropping method based on space difference.
            if view == 'sagittal':
                croprow_s = row_s - 2*buf
                croprow_f = row_f + 2*buf

            if view == 'axial':
                spacediff = 2*buf+(cropcol_f-cropcol_s)-aspect*(row_f-row_s)
                croprow_s = row_s - int(0.5*spacediff)
                croprow_f = row_f + int(0.5*spacediff)

        #if accidentally get voxel coordinate < 0, set it = 0 (don't worry about asymmetry around ROI box for now)
        if croprow_s < 0:
            croprow_s = 0
        if cropcol_s < 0:
            cropcol_s = 0

        img_wroi_crop = img_wroi[croprow_s:croprow_f, cropcol_s:cropcol_f]
        return img_wroi_crop

    buf = 5 #number of pixels of buffer around each side of larger dimension of ROI box

    #Edit 8/13/2020: Call this function to adjust values in 8-bit image using window and level values
    def adjustImgScaleWithWL(img,scale,window,level):
        #max value for window level settings is given by the equation below
        windowmax = float(level) + float(window/2)
        print("windowmax to 8 bit")
        print(str(scale*windowmax))

        #if windowmax scaled to 8 bit is less than 255, set values above this to 255
        if(scale*windowmax < 255):
            img = np.where(img < (scale*windowmax), img, 255)

        #if windowmax scaled to 8 bit is greater than 255, rescale the image so that scale*windowmax is 255
        if(scale*windowmax > 255):
            scale2 = 255/(scale*windowmax)
            img = scale2*img

        return img
            

#-------------------AXIAL---------------------------------------------------------------------------------------------------------------------------------
    #Edit 5/21/2020: Don't use subtraction images for 2nd and 3rd columns of report.
    #Instead, use early post-contrast image


    z_maxA = create2DimgAllFunctions.chooseMaxTumorSlice(zs,zf,tumor_mask) #find z-slice of axial tumor mask with max tumor area
    print("Found z slice to use, using slice #" + str(z_maxA))
    #image axial slice
    diffimg3d = img3d-preimg3d
    img_ax_slc = img3d[:,:,z_maxA] #take z slice of ROI with max tumor area

    #Edit 8/12/2020: Scale all images to max value of 255
    ax_slc_scale = 255/(np.amax(img_ax_slc))
    img_ax_slc = img_ax_slc*ax_slc_scale

    #call function to adjust image values using window level settings
    img_ax_slc = adjustImgScaleWithWL(img_ax_slc,ax_slc_scale,window,level)
    
    img_ax_slc_wroi = create2DimgAllFunctions.createImgWithROIRect(img_ax_slc,xs,xf,ys,yf,omitCount,omitradii,omitcenters,'ax',z_maxA) #np.amax(img_ax_slc)) #draw ROI rectangle on axial slice
    print("Drew rectangle on axial slice")
    print(img_ax_slc_wroi.shape)
    print("Min/Max of Axial slice array")
    print(np.amin(img_ax_slc_wroi))
    print(np.amax(img_ax_slc_wroi))

    #MIP axial slice
    mip_axial = create2DimgAllFunctions.makeMIP(diffimg3d) #make axial mip
    
    #Edit 8/12/2020: Scale all images to max value of 255
    ax_mip_scale = 255/(np.amax(mip_axial))
    mip_axial = mip_axial*ax_mip_scale

    #call function to adjust image values using window level settings
    mip_axial = adjustImgScaleWithWL(mip_axial,ax_mip_scale,window,level)
    
    mip_axial_wroi = create2DimgAllFunctions.createImgWithROIRect(mip_axial,xs,xf,ys,yf,omitCount,omitradii,omitcenters,'ax',z_maxA) #draw ROI rectangle on axial mip
    print("Made axial MIP with ROI and omit boxes")
    print(mip_axial_wroi.shape)
    print("Min/Max of Axial MIP array")
    print(np.amin(mip_axial_wroi))
    print(np.amax(mip_axial_wroi))

    #image axial slice crop with SER colormap overlayed
    #first, add SER colorization to whole axial slice
    img_ax_clr = serColorize(img_ax_slc_wroi,tumor_mask[:,:,z_maxA],ser[:,:,z_maxA],ys,yf,xs,xf)
    #then, crop colorized axial slice to small square region around ROI
    img_ax_clr_crop = sqCropROIimg(img_ax_clr,ys,yf,xs,xf,buf,1,'axial')
    print("Made axial SER colormap")
    print(img_ax_clr_crop.shape)

#-------------------SAGITTAL---------------------------------------------------------------------------------------------------------------------------------
    #Edit 5/21/2020: Don't use subtraction images for 2nd and 3rd columns of report.
    #Instead, use early post-contrast image

    #Edit 7/21/2020: Use gunzipped as source path for images if necessary
    if(gzipped == 1):
        imagepath = path + "\\" + "gunzipped"
    else:
        imagepath = path

    #5/21/2021: Using pre slice for voxels to cm^3 conversion factor to match other code

    pre_path = imagepath + "\\" + str(dce_folders[0])
    #Then, use DICOM header to compute voxel volume in cubic centimeters (cc)
    pre_imgs = os.listdir(pre_path)
    pre_imgs = sorted(pre_imgs)
    pre_img1path = pre_path + "\\" + pre_imgs[0]

    try:
      pre_hdr1 = pydicom.dcmread(pre_img1path,stop_before_pixels = True)
    except:
      pre_hdr1 = dicom.read_file(pre_img1path)

        
    #edit 6/11/2020: Split by single vs multi folder DCE instead of non-Philips vs Philips
    if (len(dce_folders) == 1):
        earlyslice1path = imagepath + "\\" + str(dce_folders[0]) + "\\" + fsort[earlyPostContrastNum][0]

    if(len(dce_folders) == 2):
        #Edit 1/29/2021: For 2 DCE folders, use earlyPostContrastNum-1 because there is no pre-contrast image in multivolume folder
        earlyslice1path = imagepath + "\\" + str(dce_folders[1]) + "\\" + fsort[earlyPostContrastNum-1][0]
        
        
    if(len(dce_folders) > 2):
        #7/6/2021: Add exception for exams where DCE series folders have number
        #and letters in name, like Duke TCIA.
        try:
            earlypath = imagepath + "\\" + str(int(dce_folders[earlyPostContrastNum]))
        except:
            earlypath = imagepath + "\\" + str(dce_folders[earlyPostContrastNum])
            
        files = os.listdir(earlypath)
        files = sorted(files)
        earlyslice1path = earlypath + "\\" + files[0]

    try:
        earlyslice1hdr = pydicom.dcmread(earlyslice1path,stop_before_pixels = True)
    except:
        earlyslice1hdr = dicom.read_file(earlyslice1path)

    #aspect ratio for sagittal images
    #Edit 5/11/2020: Using aff_mat[2,2] instead of hdr.SpacingBetweenSlices to allow aspect ratio definition for Siemens scanner images
    asp = float(abs(aff_mat[2,2]))/(float(earlyslice1hdr.PixelSpacing[1]))
    print("aspect ratio: " + str(asp))

    tumor_mask_sagittal = np.transpose(tumor_mask,(2,0,1)) #z,y,x
    ser_sagittal = np.transpose(ser,(2,0,1)) #z,y,x 
    x_maxA = create2DimgAllFunctions.chooseMaxTumorSlice(xs,xf,tumor_mask_sagittal) #find x-slice of sagittal tumor mask with max tumor area

    print("Found x slice to use, using slice #"+str(x_maxA))

    diffimg3d_sagittal = np.transpose(diffimg3d,(2,0,1)) #z,y,x
    img3d_sagittal = np.transpose(img3d,(2,0,1)) #z,y,x
    img_sag_slc = img3d_sagittal[:,:,x_maxA] #take x slice of ROI with max tumor area

    #Edit 8/12/2020: Scale all images to max value of 255
    sag_slc_scale = 255/(np.amax(img_sag_slc))
    img_sag_slc = img_sag_slc*sag_slc_scale

    #call function to adjust image values using window level settings
    img_sag_slc = adjustImgScaleWithWL(img_sag_slc,sag_slc_scale,window,level)
    
    img_sag_slc_wroi = create2DimgAllFunctions.createImgWithROIRect(img_sag_slc,ys,yf,zs,zf,omitCount,omitradii,omitcenters,'sag',x_maxA)#np.amax(img_sag_slc)) #draw ROI rectangle on sagittal slice
    print("Drew rectangle on sagittal slice")
    print(img_sag_slc_wroi.shape)
    print("Min/Max of Sagittal slice array")
    print(np.amin(img_sag_slc_wroi))
    print(np.amax(img_sag_slc_wroi))


    #Edit 6/12/2020: If x_maxA is less than halfway through x-range, use first half of x-range only for creating MIP
    #otherwise, use second half of x-range only for creating MIP
    if (x_maxA < diffimg3d_sagittal.shape[2]/2):
        diffimg3d_sagittal_formip = diffimg3d_sagittal[:,:,0:int(diffimg3d_sagittal.shape[2]/2)]
    else:
        diffimg3d_sagittal_formip = diffimg3d_sagittal[:,:,int(diffimg3d_sagittal.shape[2]/2):int(diffimg3d_sagittal.shape[2])]
        
    mip_sagittal = create2DimgAllFunctions.makeMIP(diffimg3d_sagittal_formip) #make sagittal mip

    #Edit 8/12/2020: Scale all images to max value of 255
    sag_mip_scale = 255/(np.amax(mip_sagittal))
    mip_sagittal = mip_sagittal*sag_mip_scale

    #call function to adjust image values using window level settings
    mip_sagittal = adjustImgScaleWithWL(mip_sagittal,sag_mip_scale,window,level)

    mip_sagittal_wroi = create2DimgAllFunctions.createImgWithROIRect(mip_sagittal,ys,yf,zs,zf,omitCount,omitradii,omitcenters,'sag',x_maxA) #draw ROI rectangle on sagittal mip
    print("Made sagittal MIP with roi and omit boxes")
    print(mip_sagittal_wroi.shape)
    print("Min/Max of Sagittal MIP array")
    print(np.amin(mip_sagittal_wroi))
    print(np.amax(mip_sagittal_wroi))

    #image sagittal slice crop with SER colormap overlayed
    #first, add SER colorization to whole sagittal slice
    img_sag_clr = serColorize(img_sag_slc_wroi,tumor_mask_sagittal[:,:,x_maxA],ser_sagittal[:,:,x_maxA],zs,zf,ys,yf)
    #then, crop colorized sagittal slice to small square region around ROI
    img_sag_clr_crop = sqCropROIimg(img_sag_clr,zs,zf,ys,yf,buf,asp,'sagittal')
    #flip l/r to give desired orientation
    img_sag_clr_crop = np.fliplr(img_sag_clr_crop)
    print("Made sagittal SER colormap")
    print(img_sag_clr_crop.shape)

    #Now that all outputs have been generated, crop extra space off of sagittal slice and MIP to make less oblong
    #If done earlier, may get error while creating sagittal SER colorization
    img_sag_slc_wroi = np.fliplr(img_sag_slc_wroi)
    mip_sagittal_wroi = np.fliplr(mip_sagittal_wroi)

#------------------------------creating figure------------------------------------------------------------------------------------------------------------

    #report JPEG's name is save as report PDF's name but with extension changed
    savenamejpeg = savenamepdf[:-3] + "jpeg"
    savenamedcm = savenamepdf[:-3] + "dcm"

    #Edit 5/13/2020: Always want the image orientations seen in the GE reports. So if other exams are not already like that,
    #apply fliplr or flipud so that you can keep the orientation label positions.
    
    #Edit 4/30/2020: Use affine matrix values to determine appropriate orientation labels on images
    
    if (float(aff_mat[0,1])>0):
        mip_axial_wroi = np.fliplr(mip_axial_wroi)
        img_ax_slc_wroi = np.fliplr(img_ax_slc_wroi)
        img_ax_clr_crop = np.fliplr(img_ax_clr_crop)
       
    axleftlbl = 'L'
    axrightlbl = 'R'

    #for sagittal, A/P labels are reversed from what affine matrix implies because we have applied
    #a lr flip to the sagittal images
    if (float(aff_mat[1,0])>0):
        mip_axial_wroi = np.flipud(mip_axial_wroi)
        img_ax_slc_wroi = np.flipud(img_ax_slc_wroi)
        img_ax_clr_crop = np.flipud(img_ax_clr_crop)

        #for these images, this would be the 2nd fliplr, back to original l/r orientation
        mip_sagittal_wroi = np.fliplr(mip_sagittal_wroi)
        img_sag_slc_wroi = np.fliplr(img_sag_slc_wroi)
        img_sag_clr_crop = np.fliplr(img_sag_clr_crop)
        
    axtoplbl = 'P'
    axbtmlbl = 'A'

    sagleftlbl = 'A'
    sagrightlbl = 'P'



    #Edit 6/26/2020: Using original orientation from:
    #loadables = plugin.examine([dicom_names])
    #volume = plugin.load(loadables[0])
    #Since this reads slices in reverse order from what is expected in aff_mat, the logic below
    #must be inverted
    if (float(aff_mat[2,2])<0):
        mip_sagittal_wroi = np.flipud(mip_sagittal_wroi)
        img_sag_slc_wroi = np.flipud(img_sag_slc_wroi)
        img_sag_clr_crop = np.flipud(img_sag_clr_crop)
        
#Edit 6/30/2020: No more imgReversed
    #If you had to invert axial slices of ROI, img slices are reversed from what aff_mat suggests and you have
    #to compensate for this with a flip
##    if(imgReversed == 1):
##        mip_sagittal_wroi = np.flipud(mip_sagittal_wroi)
##        img_sag_slc_wroi = np.flipud(img_sag_slc_wroi)
##        img_sag_clr_crop = np.flipud(img_sag_clr_crop)


    #Edit 6/29/2020: Compare z value of intial position in IJKToRAS matrix from inputVolume to same element in
    #my IJKToLPS aff_mat, and if they do match, add an extra flip to sagittal images in report
    #If these two z origin values match, this means Slicer used closest to foot as 1st slice and
    #closest to head as last slice, which is why you have to flip
    zorig = ijkToRASmat.GetElement(2,3)
    zorig_LPS = aff_mat[2,3]
    if round(zorig,2) == round(zorig_LPS,2):
        mip_sagittal_wroi = np.flipud(mip_sagittal_wroi)
        img_sag_slc_wroi = np.flipud(img_sag_slc_wroi)
        img_sag_clr_crop = np.flipud(img_sag_clr_crop)
        

    sagtoplbl = 'H'
    sagbtmlbl = 'F'



    #read header of slice in late post-contrast image
    #edit 6/11/2020: Split by single vs multi folder DCE instead of non-Philips vs Philips
    if (len(dce_folders) == 1):
        lateslice1path = imagepath + "\\" + str(dce_folders[0]) + "\\" + fsort[latePostContrastNum][0]

    if(len(dce_folders) == 2):
        #Edit 1/29/2021: For 2 DCE folders, use latePostContrastNum-1 because there is no pre-contrast image in multivolume folder
        lateslice1path = imagepath + "\\" + str(dce_folders[1]) + "\\" + fsort[latePostContrastNum-1][0]
        
    if(len(dce_folders) > 2):
        #7/6/2021: Add exception for exams where DCE series folders have number
        #and letters in name, like Duke TCIA.
        try:
            latepath = imagepath + "\\" + str(int(dce_folders[latePostContrastNum]))
        except:
            latepath = imagepath + "\\" + str(dce_folders[latePostContrastNum])
        files = os.listdir(latepath)
        files = sorted(files)
        lateslice1path = latepath + "\\" + files[0]

    try:
        lateslice1hdr = pydicom.dcmread(lateslice1path,stop_before_pixels = True)
    except:
        lateslice1hdr = dicom.read_file(lateslice1path)


    font = {'color': 'yellow','size': 15} #font dictionary settings for orientation labels on images
    fontfov = {'color': 'yellow', 'size': 8} #font dictionary settings for FOV labels on images
    fontfovcrop = {'color': 'yellow', 'size': 6} #font dictionary settings for FOV labels on images

#--------------Create strings that will go in figure---------------------------#
    patient_id = "Patient ID: " + idstr #Using Patient ID (I-SPY ID)

    #7/19/2021: Create string for manufacturer
    manufact_str = "Manufacturer: " + manufacturer

    try:
        institution = "Institution: " + str(earlyslice1hdr.ClinicalTrialSiteName) #Using Clinical Trial Site Name
    except: #in some cases, have to use Institution Name field instead
        try:
            institution = "Institution: " + str(earlyslice1hdr[0x8,0x80].value)
        except:
            institution = 'Site Unknown'
            
    #preparing study string using early post-contrast image header
    studytime = str(earlyslice1hdr.StudyTime)
    studydate = str(earlyslice1hdr.StudyDate)
    study = "Study: "+studydate[4:6]+"/"+studydate[6:]+"/"+studydate[0:4]+" "+studytime[0:2]+":"+studytime[2:4]+":"+studytime[4:]

    #visit id
    #Edit 7/30/2020: Use nodevisstr, which was obtained from exampath instead of DICOM header, to fill this field
    visit_id = "Visit ID: " + nodevisstr
##    try:
##        visit_id = "Visit ID: " + str(earlyslice1hdr.ClinicalTrialTimePointID)
##    except: #in some cases, this is not in the DICOM header
##        visit_id = "Visit ID: Unavailable"

    #Return Left or Right breast based on xcenter of ROI
    #Edit 5/1/2020: make it so that L/R determination doesn't depend
    #on previously existing xml file
    xcen = (xs+xf)/2
    imgcen = img_sag_slc_wroi.shape[1]/2
    if (xcen < imgcen and float(aff_mat[0,1])<0) or (xcen > imgcen and float(aff_mat[0,1])>0):
        breast = "Breast: L"
    else:
        breast = "Breast: R"

    
##    centercoords, widthcoords, heightcoords, depthcoords = parse_xml.returnPrimaryROICoordsFromXML(path)
##    if (int(centercoords[0])>0):
##        breast = "Breast: L"
##    else:
##        breast = "Breast: R"

    #tumor volume
    #Counting all white voxels in VOI as part of tumor
    tumor_mask_roi = tumor_mask*voi_mask
    #tumor_mask_roi = tumor_mask_roi[ys:yf,xs:xf,zs:zf]
    ftv_voxels = np.sum(tumor_mask_roi)
    #Edit 5/11/2020: Using aff_mat[2,2] instead of hdr.SpacingBetweenSlices to allow aspect ratio definition for Siemens scanner images
    voxsize_mm3 = float(pre_hdr1.PixelSpacing[0])*float(pre_hdr1.PixelSpacing[1])*float(abs(aff_mat[2,2]))
    voxsize_cm3 = voxsize_mm3/1000
    ftv_cm3 = ftv_voxels*voxsize_cm3
    tumor_volume = "Tumor Volume: {} cc".format(round(ftv_cm3,3))
    
    #auto timing
    #Format:
    # MM:SS from content time of early (early #) / MM:SS from content time of late (late #)
    earlycontenttime = str(earlyslice1hdr.ContentTime)
    latecontenttime = str(lateslice1hdr.ContentTime)
    
    if earlydiffss<10:
        earlydiffss_str = "0"+str(earlydiffss)
    else:
        earlydiffss_str = str(earlydiffss)

    if latediffss<10:
        latediffss_str = "0"+str(latediffss)
    else:
        latediffss_str = str(latediffss)

    autotiming = "Auto Timing: " + str(earlydiffmm) + ":" + earlydiffss_str + "(" + str(earlyPostContrastNum) + ")" + "/" + str(latediffmm) + ":" + latediffss_str + "(" + str(latePostContrastNum) + ")"

    #roibounds
    roi = "ROI: [{},{},{}]:[{},{},{}]".format(xs,ys,zs,xf,yf,zf)
    #pe and min connected neighbors thresholds
    pe_mnc = "PE/MNC Threshold: {}/{}".format(pethresh,minconnpix)

    #Edit 6/12/2020: Use tempres from chooseEarlyLate code, which is already in seconds
    scanduration = "Scan Duration: " + str(round(float(tempres),1)) + " s"
    
    #pre-contrast threshold (voxel value and % of max (or 95%ile value) in ROI 
    gray_thresh = "Gray Threshold: {}/{}%".format(round(pre_thresh,1),100*pct)

    #create strings that show FOV for each image
    ax_fovx = round(earlyslice1hdr.PixelSpacing[0]*img_ax_slc_wroi.shape[1]/10,2)
    ax_fovy = round(earlyslice1hdr.PixelSpacing[1]*img_ax_slc_wroi.shape[0]/10,2)
    ax_full_fov = "FOV: " + str(ax_fovx) + " x " + str(ax_fovy) + " cm"

    ax_crop_fovx = round(earlyslice1hdr.PixelSpacing[0]*img_ax_clr_crop.shape[1]/10,2)
    ax_crop_fovy = round(earlyslice1hdr.PixelSpacing[1]*img_ax_clr_crop.shape[0]/10,2)
    ax_crop_fov = "FOV: " + str(ax_crop_fovx) + " x " + str(ax_crop_fovy) + " cm"

    sag_fovx = round(earlyslice1hdr.PixelSpacing[1]*img_sag_slc_wroi.shape[1]/10,2)
    sag_fovy = abs(round(aff_mat[2,2]*img_sag_slc_wroi.shape[0]/10,2))
    sag_full_fov = "FOV: " + str(sag_fovx) + " x " + str(sag_fovy) + " cm"

    sag_crop_fovx = round(earlyslice1hdr.PixelSpacing[1]*img_sag_clr_crop.shape[1]/10,2)
    sag_crop_fovy = abs(round(aff_mat[2,2]*img_sag_clr_crop.shape[0]/10,2))
    sag_crop_fov = "FOV: " + str(sag_crop_fovx) + " x " + str(sag_crop_fovy) + " cm"

    
    #Read UCSF logo into array
    #Edit 7/30/2020: Use 3D Slicer logo instead because Jessica said use of UCSF logo
    #is confusing for report of exam from non-UCSF site
    #7/8/2021: Make this part work no matter which directory & computer
    #the modules are stored in.
    #The full path to ftv_plots.py is
    #automatically stored in __file__ so I'm using that.
    pathtofunc,funcname = os.path.split(os.path.realpath(__file__))
    logo_path = os.path.join(pathtofunc,'3DSlicerLogo.png')
    logo = cv2.imread(logo_path)


    
#------------------------Plotting starts here----------------------------------#
    fig = plt.figure(constrained_layout=True,figsize=(11,8.5))
    spec = gridspec.GridSpec(ncols=3,nrows=5,height_ratios=[1,3,4,1,0.5],figure=fig) #this height ratio is chosen because axial image is ~3x taller than sagittal, but still want sagittal to look big
                                                                    #width_ratios=[3,3,1],This width ratio is chosen because VOI width is between 1/5 and 1/3 of regular 512 pixel width, but want VOI zoom to look big 
    #Top text content
    #top left
    a = fig.add_subplot(spec[0,0])
    txtplot = plt.text(0,0.6,patient_id,fontsize = 10)
    txtplot = plt.text(0,0.2,study,fontsize=10)
    plt.axis('off')
    #top center
    a = fig.add_subplot(spec[0,1])
    txtplot = plt.text(0,0.7,visit_id,fontsize = 10)
    txtplot = plt.text(0,0.4,institution,fontsize = 10)
    txtplot = plt.text(0,0.1,manufact_str,fontsize = 10)
    plt.axis('off')
    #top right
    a = fig.add_subplot(spec[0,2])
    imgplot = plt.imshow(logo)
    plt.axis('off')

    #Edit 4/29/2020: Finalized orientation label letters on images
    
    #Sagittal subtraction MIP
    a = fig.add_subplot(spec[1,0])
    imgplot = plt.imshow(mip_sagittal_wroi,cmap='gray',aspect=asp)#,vmin=dispmin/sag_mip_normfact,vmax=dispmax/sag_mip_normfact)
    a.set_title('Sagittal - Subtracted MIP')
    txtfov = plt.text( (3*mip_sagittal_wroi.shape[1]/5), (49*mip_sagittal_wroi.shape[0]/50), sag_full_fov, fontdict=fontfov) #FOV in bottom right
    plt.axis('off')
    
    #Sagittal regular MRI center x-slice of ROI
    a = fig.add_subplot(spec[1,1])
    imgplot = plt.imshow(img_sag_slc_wroi,cmap='gray',aspect=asp)#,vmin=dispmin/img_sag_slc_normfact,vmax=dispmax/img_sag_slc_normfact)
    txtplot = plt.text((img_sag_slc_wroi.shape[1]/50),(img_sag_slc_wroi.shape[0]/2),sagleftlbl,fontdict=font) # 1/50 through column range, 1/2 through row range of image
    txtplot = plt.text((47*img_sag_slc_wroi.shape[1]/50),(img_sag_slc_wroi.shape[0]/2),sagrightlbl,fontdict=font) # 47/50 through column range, 1/2 through row range of image
    txtplot = plt.text((img_sag_slc_wroi.shape[1]/2),(5*img_sag_slc_wroi.shape[0]/50),sagtoplbl,fontdict=font) #1/2 through column range, 5/50 through row range of image
    txtplot = plt.text((img_sag_slc_wroi.shape[1]/2),(49*img_sag_slc_wroi.shape[0]/50),sagbtmlbl,fontdict=font) #1/2 through column range, 49/50 through row range of image
    txtfov = plt.text( (3*img_sag_slc_wroi.shape[1]/5), (49*img_sag_slc_wroi.shape[0]/50), sag_full_fov, fontdict=fontfov) #FOV in bottom right
    a.set_title('Sagittal')
    plt.axis('off')
    
    #Sagittal cropped slice with SER colorization
    a = fig.add_subplot(spec[1,2])
    imgplot = plt.imshow(img_sag_clr_crop,aspect=asp)
    txtplot = plt.text((img_sag_clr_crop.shape[1]/50),(img_sag_clr_crop.shape[0]/2),sagleftlbl,fontdict=font) # 1/50 through column range, 1/2 through row range of image
    txtplot = plt.text((47*img_sag_clr_crop.shape[1]/50),(img_sag_clr_crop.shape[0]/2),sagrightlbl,fontdict=font) # 47/50 through column range, 1/2 through row range of image
    txtplot = plt.text((img_sag_clr_crop.shape[1]/2),(5*img_sag_clr_crop.shape[0]/50),sagtoplbl,fontdict=font) #1/2 through column range, 5/50 through row range of image
    txtplot = plt.text((img_sag_clr_crop.shape[1]/2),(49*img_sag_clr_crop.shape[0]/50),sagbtmlbl,fontdict=font) #1/2 through column range, 49/50 through row range of image
    txtfov = plt.text( (3*img_sag_clr_crop.shape[1]/5), (49*img_sag_clr_crop.shape[0]/50), sag_crop_fov, fontdict=fontfovcrop) #FOV in bottom right
    a.set_title('Sagittal - Lesion Focus')
    plt.axis('off')
        
    #Axial subtraction MIP
    a = fig.add_subplot(spec[2,0])
    imgplot = plt.imshow(mip_axial_wroi,cmap='gray')#,vmin=dispmin/ax_mip_normfact,vmax=dispmax/ax_mip_normfact)
    a.set_title('Axial - Subtracted MIP')
    txtfov = plt.text( (3*mip_axial_wroi.shape[1]/5), (49*mip_axial_wroi.shape[0]/50), ax_full_fov, fontdict=fontfov) #FOV in bottom right
    plt.axis('off')
    #Axial regular MRI center z-slice of ROI
    a = fig.add_subplot(spec[2,1])
    imgplot = plt.imshow(img_ax_slc_wroi,cmap='gray')#,vmin=dispmin/img_ax_slc_normfact,vmax=dispmax/img_ax_slc_normfact)
    txtplot = plt.text((img_ax_slc_wroi.shape[1]/50),(img_ax_slc_wroi.shape[0]/2),axleftlbl,fontdict=font) #1/50 through column range, 1/2 through row range of image
    txtplot = plt.text((47*img_ax_slc_wroi.shape[1]/50),(img_ax_slc_wroi.shape[0]/2),axrightlbl,fontdict=font) #47/50 through column range, 1/2 through row range of image
    txtplot = plt.text((img_ax_slc_wroi.shape[1]/2),(5*img_ax_slc_wroi.shape[0]/50),axtoplbl,fontdict=font) #1/2 through column range, 5/50 through row range of image
    txtplot = plt.text((img_ax_slc_wroi.shape[1]/2),(49*img_ax_slc_wroi.shape[0]/50),axbtmlbl,fontdict=font) #1/2 through column range, 49/50 through row range of image
    txtfov = plt.text( (3*img_ax_slc_wroi.shape[1]/5), (49*img_ax_slc_wroi.shape[0]/50), ax_full_fov, fontdict=fontfov) #FOV in bottom right
    a.set_title('Axial')
    plt.axis('off')
    #Axial cropped slice with SER colorization
    a = fig.add_subplot(spec[2,2])
    imgplot = plt.imshow(img_ax_clr_crop)
    txtplot = plt.text((img_ax_clr_crop.shape[1]/50),(img_ax_clr_crop.shape[0]/2),axleftlbl,fontdict=font) #1/50 through column range, 1/2 through row range of image
    txtplot = plt.text((47*img_ax_clr_crop.shape[1]/50),(img_ax_clr_crop.shape[0]/2),axrightlbl,fontdict=font) #47/50 through column range, 1/2 through row range of image
    txtplot = plt.text((img_ax_clr_crop.shape[1]/2),(5*img_ax_clr_crop.shape[0]/50),axtoplbl,fontdict=font) #1/2 through column range, 5/50 through row range of image
    txtplot = plt.text((img_ax_clr_crop.shape[1]/2),(49*img_ax_clr_crop.shape[0]/50),axbtmlbl,fontdict=font) #1/2 through column range, 49/50 through row range of image
    txtfov = plt.text( (3*img_ax_clr_crop.shape[1]/5), (49*img_ax_clr_crop.shape[0]/50), ax_crop_fov, fontdict=fontfovcrop) #FOV in bottom right
    a.set_title('Axial - Lesion Focus')
    plt.axis('off')
    
    #Text region 1
    a = fig.add_subplot(spec[3,0])
    a.set_title('Findings')
    txtplot = plt.text(0.1,0.9,breast,fontsize=12)
    txtplot = plt.text(0.1,0.5,tumor_volume,fontsize=12)
    plt.axis('off')
    #Text region 2
    a = fig.add_subplot(spec[3,1])
    txtplot = plt.text(0.1,1,'Settings',fontsize = 15)
    txtplot = plt.text(0.1,0.8,autotiming,fontsize=11)
    txtplot = plt.text(0.1,0.6,pe_mnc,fontsize=11)
    txtplot = plt.text(0.1,0.4,scanduration,fontsize=11)
    txtplot = plt.text(0.1,0.2,gray_thresh,fontsize=11)
    txtplot = plt.text(0.1,0,roi,fontsize=11)
    plt.axis('off')
    #Text region 3
    a = fig.add_subplot(spec[3,2])
    txtplot = plt.text(0.1,0.9,'__________________________________',fontsize = 11)
    txtplot = plt.text(0.1,0.65,'Radiologist Signature',fontsize = 12)
    txtplot = plt.text(0.1,0.25,'__________________________________',fontsize = 11)
    txtplot = plt.text(0.1,0,'Date',fontsize = 12)
    plt.axis('off')

    #footer info
    #bottom left
    a = fig.add_subplot(spec[4,0])
    txtplot = plt.text(0.1,0.7,'NOT FOR DIAGNOSIS', fontsize = 15)
    plt.axis('off')
    #adding approval info to bottom center
    a = fig.add_subplot(spec[4,1])
    txtplot = plt.text(0,1,'Approval',fontsize = 15)
    txtplot = plt.text(0.06,0.5,'Original Image Reviewed BI-RADS Readable',fontsize = 11)
    txtplot = plt.text(0.06,0.1,'Volume Verified and Acceptable',fontsize = 11)
    rect1 = mpatches.Rectangle((0,0.5),width=0.025,height=0.2,fill=False) #lower left at (0.1,0.5), width & height 0.15
    rect2 = mpatches.Rectangle((0,0.1),width=0.025,height=0.2,fill=False) #lower left at (0.1,0.1), width & height 0.15
    a.add_patch(rect1)
    a.add_patch(rect2)
    plt.axis('off')
    #bottom right
    a = fig.add_subplot(spec[4,2])
    txtplot = plt.text(0.5,0.7,'Page 1 of 1',fontsize = 11)
    plt.axis('off')
    
    #plt.tight_layout()
    fig.savefig(savenamejpeg,orientation='landscape', dpi = 300, bbox_inches=0) #save as jpeg
    fig.savefig(savenamepdf,orientation='landscape', dpi = 300, bbox_inches=0) #save as pdf

    #save as dicom too
    command = "img2dcm" + " " + savenamejpeg + " " + savenamedcm #use img2dcm to create dicom report from jpeg
    os.system(command)

    return
