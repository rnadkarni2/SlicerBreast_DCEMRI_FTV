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
#This code contains the function that creates
#the SER colorized lesion segmentation after the user
#clicks on the "Show Lesion Segmentation" checkbox
#in module 2.


#Function to colorize tumor in ROI based on SER values
def serColorize3D(rgbimg,tumor_mask,voi_mask,ser,zs,zf,ys,yf,xs,xf):
    #only interested in colorization where tumor mask and voi mask have value 1
    tumor_voi_mask = tumor_mask*voi_mask
    nblue = 0
    npurple = 0
    ngreen = 0
    nred = 0
    nyellow = 0

    #Only doing colorization within ROI
    for z in range(zs,zf+1):
        for y in range(ys,yf+1):
            for x in range(xs,xf+1):
                #Only add color to voxels that have been labeled as tumor
                if tumor_voi_mask[x,y,z] == 1:
                    #If SER in range 0 to 0.9, make the voxel blue
                    if ser[x,y,z] > 0 and ser[x,y,z] <= 0.9:
                        rgbimg[x,y,z,:] = [0,0,255]
                        nblue = nblue+1 #update # of blue voxels

                    #If SER in range 0.9 to 1, make the voxel purple
                    if ser[x,y,z] > 0.9 and ser[x,y,z] <= 1:
                        rgbimg[x,y,z,:] = [128,0,128]
                        npurple = npurple+1 #update # of purple voxels

                    #If SER in range 1 to 1.3, make the voxel green
                    if ser[x,y,z] > 1 and ser[x,y,z] <= 1.3:
                        rgbimg[x,y,z,:] = [0,255,0]
                        ngreen = ngreen+1 #update # of green voxels

                    #If SER in range 1.3 to 1.75, make the voxel red
                    if ser[x,y,z] > 1.3 and ser[x,y,z] <= 1.75:
                        rgbimg[x,y,z,:] = [255,0,0]
                        nred = nred+1 #update # of red voxels

                    #If SER above 1.75, make the voxel yellow
                    #Edit 5/1/2020: Make last one yellow instead of white to avoid
                    #confusion with bright pixels in grayscale image and cap at max
                    #value 3 as seen in Aegis
                    if ser[x,y,z] > 1.75 and ser[x,y,z] <= 3:
                        rgbimg[x,y,z,:] = [255,255,0]
                        nyellow = nyellow+1 #update # of yellow voxels


                    #Edit 7/10/2020: Now that I am able to overlay this on top of other image,
                    #make all non SER colored voxels white, instead of same as background image
                    if ser[x,y,z] <= 0 or ser[x,y,z] > 3:
                        rgbimg[x,y,z,:] = [255,255,255]

    #Convert # of tumor of voxels of each color into % of total # of tumor voxels
    pctblue = 100*float(nblue)/(float(nblue+npurple+ngreen+nred+nyellow))
    pctpurple = 100*float(npurple)/(float(nblue+npurple+ngreen+nred+nyellow))
    pctgreen = 100*float(ngreen)/(float(nblue+npurple+ngreen+nred+nyellow))
    pctred = 100*float(nred)/(float(nblue+npurple+ngreen+nred+nyellow))
    pctyellow = 100*float(nyellow)/(float(nblue+npurple+ngreen+nred+nyellow))

    return rgbimg, nblue, pctblue, npurple, pctpurple, ngreen, pctgreen, nred, pctred, nyellow, pctyellow
