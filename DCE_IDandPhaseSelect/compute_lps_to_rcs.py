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



#Rohan Nadkarni
#File for reading the xml file specifying VOI and omit regions in LPS coordinates and outputting voxel range for these regions

import os
import pydicom
import dicom
import numpy as np
import re

def computeAffineAndAffineInverse(exampath,prefoldernum,nslice,fsort):
    imgpath = exampath+"\\"+str(prefoldernum)
    print(imgpath)
    files = os.listdir(imgpath)

    #nslice = 0 is when there is 1 phase per folder, so N = slices/phase = # of DICOMs in folder
    #If all DCE in same folder, nslice is already set to slices/phase, so N = nslice
    if (nslice == 0):
        N = len(files)
        
        file1search1 = [i for i in files if '001.dcm' in i]
        file1search2 = [i for i in files if '001.DCM' in i]

        #Added these 2 due to file naming in UCSF ISPY ID 16078
        file1search3 = [i for i in files if 'I1.dcm' in i]
        file1search4 = [i for i in files if 'I1.DCM' in i]

        if len(file1search1) > 0:
            file1 = imgpath + "\\" + file1search1[0]

        if len(file1search2) > 0:
            file1 = imgpath + "\\" + file1search2[0]

        if len(file1search3) > 0:
            file1 = imgpath + "\\" + file1search3[0]
            
        if len(file1search4) > 0:
            file1 = imgpath + "\\" + file1search4[0]

        file1search5 = []
        #Edit 1/26/2021: file1search that incorporates DICOMs
        #with no .DCM or .dcm extension
        if( len(file1search1)==0 and len(file1search2)==0 and len(file1search3)==0 and len(file1search4)==0):
            file1search5 = [i for i in files if i.isdigit()]

            if len(file1search5) > 0:
                file1search5 = sorted(file1search5)
                file1 = imgpath + "\\" + file1search5[0]
            
    else:
        N = nslice
        file1 = imgpath + "\\" + fsort[0][0]


 
    try:
        img1 = pydicom.dcmread(file1)
    except:
        img1 = dicom.read_file(file1)

    img1_orient = img1[0x20,0x37] #Image Orientation (Patient)
    img1_orient = [float(i) for i in img1_orient] #convert from list of strings to list of floats (numeric)

    x1_orient = img1_orient[0:3] #x,y,z change as you move 1 column to the right

    y1_orient = img1_orient[3:] #x,y,z change as you move 1 row down

    img1_pos = img1[0x20,0x32] #Image Position (Patient) for 1st slice
    #Next 2 lines convert from header field with values to list of floats
    img1_pos = img1_pos[0:]
    img1_pos = [float(i) for i in img1_pos] #convert from list of strings to list of floats (numeric)

    img_sp = img1[0x28,0x30] #Pixel Spacing
    row_sp = float(img_sp[0]) #row spacing
    col_sp = float(img_sp[1]) #column spacing

    #Read DICOM for last slice
    #Once again, separate by Philips and non-Philips (6/9/2020)
    if (nslice == 0):
        if (len(file1search1)>0 or len(file1search3)>0):
            searchstr = str(N) + '.dcm'
        if (len(file1search2)>0 or len(file1search4)>0):
            searchstr = str(N) + '.DCM'

        if( len(file1search5) > 0):
            fileN = imgpath + "\\" + file1search5[len(file1search5) - 1]
        else:
            fileNsearch = [i for i in files if searchstr in i]
            fileN = imgpath + "\\" + fileNsearch[0]
    else:
        fileN = imgpath + "\\" + fsort[0][N-1]
    
    try:
        imgN = pydicom.dcmread(fileN)
    except:
        imgN = dicom.read_file(fileN)

    imgN_pos = imgN[0x20,0x32] #Image Position (Patient) for last slice
    #Next 2 lines convert from header field with values to list of floats
    imgN_pos = imgN_pos[0:]
    imgN_pos = [float(i) for i in imgN_pos] #convert from list of strings to list of floats (numeric)

    #Construct affine matrix for ijk to LPS transform using DICOM header information
    aff_mat = np.zeros((4,4)) #affine matrix for ijk to LPS transform
    aff_mat[0:3,0] = np.transpose(row_sp*np.array(y1_orient))
    aff_mat[0:3,1] = np.transpose(col_sp*np.array(x1_orient))
    aff_mat[0:3,2] = np.transpose((np.array(img1_pos)-np.array(imgN_pos))/(1-N))
    aff_mat[0:3,3] = np.transpose(np.array(img1_pos))
    aff_mat[3,3] = 1

    aff_inv_mat = np.linalg.inv(aff_mat) #inverse of affine matrix for LPS to ijk transform

    return aff_mat, aff_inv_mat




def getVOIVoxelsFromInverseAffine(exampath,xmlfilepath,prefoldernum,nslice,fsort):
    print("-----Starting conversion of VOI to voxel coordinates-----")
    aff_mat, aff_inv_mat = computeAffineAndAffineInverse(exampath,prefoldernum,nslice,fsort)

    #Edit 5/19/2020: comment out os.listdir part because you're just
    #using the xmlfilepath specified by user
    
    #Now, convert extreme ends of VOI from xml file from LPS to JIK to get voxel range of VOI
    #voipath = exampath + "\\" + "voi_lps_files"
    #xml_files = [f for f in os.listdir(voipath) if f.endswith('.xml')]
    #xml_file = xml_files[0]

    #filepath = voipath + "\\" + xml_file
    #This method successfully read all lines in the file
    with open (xmlfilepath,"r") as myfile:
        data = myfile.readlines()

    def getcoords(xmlline):
        temp = re.findall(r"[-+]?\d*\.\d+|\d+",xmlline) #find all positive or negative numbers with decimal in the line
        coords = list(map(float,temp))
        return coords

    def getXYZcoords(data,datarow):
        coordsrow = data[datarow]
        coords = getcoords(coordsrow)
        coords = np.array(coords)
        #6/9/2020: add exception for case where it splits x^-10 into x, 10
        if len(coords) > 3:
            coords_orig = coords
            coords = []
            for i in range(len(coords_orig)):
                if ( (coords_orig[i]).is_integer() == 0 or coords_orig[i] == 0): #if it's a positive integer, assume it is an exponent (this if statement is reverse of that)
                    coords.append(coords_orig[i])
                else:
                    coords[len(coords)-1] = float(coords_orig[i-1])*(10**(-int(coords_orig[i])))
            coords = np.array(coords)
        return coords

    def findWordIndices(word, data):
        wordlen = len(word)
        wordind = []
        for i in range(len(data)):
            currline = data[i]
            if (word in currline):
                wordind.append(i)
        return wordind

    #find which lines contain the relevant coordinates & check if there are >1 VOI's
    centerinds = findWordIndices('<Center',data)
    widthinds = findWordIndices('<WidthHalfLength',data)
    heightinds = findWordIndices('<HeightHalfLength',data)
    depthinds = findWordIndices('<DepthHalfLength',data)


    #First step: Do coordinate conversion from LPS to RCS for actual VOI
    centercoords = getXYZcoords(data,centerinds[0])
    print("center")
    print(centercoords)

    radius = np.zeros((3,1))
    radius.astype(float)
    widthcoords = getXYZcoords(data,widthinds[0])
    print("width")
    print(widthcoords)
    radius[0] = widthcoords[0]

    heightcoords = getXYZcoords(data,heightinds[0])
    radius[1] = heightcoords[1]
    print("height")
    print(heightcoords)

    depthcoords = getXYZcoords(data,depthinds[0])
    radius[2] = depthcoords[2]
    print("depth")
    print(depthcoords)





    def lpsTorcs(lpscoords,aff_inv_mat):
        lps_homog = np.zeros((4,1))
        lps_homog[3,0] = 1
        lps_homog[0:3,0] = np.transpose(np.array(lpscoords))
        rcs_homog = np.matmul(aff_inv_mat,lps_homog) #row,column,slice homogeneous coordinates
        return rcs_homog

    def getRCSrange(centercoords,widthcoords,heightcoords,depthcoords,aff_inv_mat):
        corner1lps = centercoords-np.abs(widthcoords+heightcoords+depthcoords)
        corner2lps = centercoords+np.abs(widthcoords+heightcoords+depthcoords)

        corner1rcs = lpsTorcs(corner1lps,aff_inv_mat)
        yf = int(corner1rcs[0])
        xf = int(corner1rcs[1])
        zf = int(corner1rcs[2])

        corner2rcs = lpsTorcs(corner2lps,aff_inv_mat)
        ys = int(corner2rcs[0])
        xs = int(corner2rcs[1])
        zs = int(corner2rcs[2])
        return xs, xf, ys, yf, zs, zf

    #Call getRCSrange for the tumor VOI
    xs, xf, ys, yf, zs, zf = getRCSrange(centercoords,widthcoords,heightcoords,depthcoords,aff_inv_mat)
    print("-----Computed VOI voxel coordinates-----")

    imgpath = exampath+"\\"+str(prefoldernum)
    files = os.listdir(imgpath)
    N = len(files)

    #Read DICOM for 1st slice
    file1 = imgpath+"\\"+files[0]
    try:
        img1 = pydicom.dcmread(file1)
    except:
        img1 = dicom.read_file(file1)

    voi_mask = np.zeros((int(img1.Rows),int(img1.Columns),len(files)))
    voi_mask[xs:xf,ys:yf,zs:zf] = 1
    #Second step: If there are omit regions, do conversion from LPS to RCS for each of those
    if(len(centerinds)>1):
        #arrays of xs, xf, ys, yf, zs, zf for all omit regions, where each index corresponds to a specific omit region
        #eg for 1st omit region, take 0th index of all 6 arrays to get its voxel range
        oxs = np.zeros((len(centerinds)-1,1))
        oxf = np.zeros((len(centerinds)-1,1))

        oys = np.zeros((len(centerinds)-1,1))
        oyf = np.zeros((len(centerinds)-1,1))

        ozs = np.zeros((len(centerinds)-1,1))
        ozf = np.zeros((len(centerinds)-1,1))
                            
        for i in range(1,len(centerinds)): #loop through all omit regions
            #center and half lengths in all directions for omit region of current iteration
            ocentercoords = getXYZcoords(data,centerinds[i])
            owidthcoords = getXYZcoords(data,widthinds[i])
            oheightcoords = getXYZcoords(data,heightinds[i])
            odepthcoords = getXYZcoords(data,depthinds[i])

            #get voxel coords of omit region and 0 out this region in VOI mask
            oxs[i-1],oxf[i-1],oys[i-1],oyf[i-1],ozs[i-1],ozf[i-1] = getRCSrange(ocentercoords,owidthcoords,oheightcoords,odepthcoords,aff_inv_mat)
            voi_mask[int(oxs[i-1]):int(oxf[i-1]),int(oys[i-1]):int(oyf[i-1]),int(ozs[i-1]):int(ozf[i-1])] = 0
    else:
        oxs = np.zeros((1,1))
        oxs[0] = -1

        oxf = np.zeros((1,1))
        oxf[0] = -1

        oys = np.zeros((1,1))
        oys[0] = -1

        oyf = np.zeros((1,1))
        oyf[0] = -1


        ozs = np.zeros((1,1))
        ozs[0] = -1

        ozf = np.zeros((1,1))
        ozf[0] = -1

    return xs, xf, ys, yf, zs, zf, oxs, oxf, oys, oyf, ozs, ozf, voi_mask, aff_mat, aff_inv_mat


#6/25/2020: Function to return RAS coords in case I decide to change coordinate transform methods to
#better handle issue of slices being read in wrong order
def readRASCoordsFromXMLFile(xmlfilepath):
    with open (xmlfilepath,"r") as myfile:
        data = myfile.readlines()

    def getcoords(xmlline):
        temp = re.findall(r"[-+]?\d*\.\d+|\d+",xmlline) #find all positive or negative numbers with decimal in the line
        coords = list(map(float,temp))
        return coords

    def getXYZcoords(data,datarow):
        coordsrow = data[datarow]
        coords = getcoords(coordsrow)
        coords = np.array(coords)

        #7/20/2020: Make compatible with xml files where (for ex) center, xcoord, ycoord, and zcoord are all on separate lines
        #I tried this for UCSD ispy2 exam and x and y range of ROI pixels didn't match with brtool, although z range did.
        #Leave this part as is--the results were fine for an Emory exam (72455) with this kind of xml file, so the problem is specific
        #to that UCSD exam (03833)
        if(len(coords) == 0):
            coords = []

            #read in x coordinate
            xrow = datarow+1
            xcoordsrow = data[xrow]
            xcoords = getcoords(xcoordsrow)
            coords.append(float(xcoords[0]))

            #read in y coordinate
            yrow = datarow+2
            ycoordsrow = data[yrow]
            ycoords = getcoords(ycoordsrow)
            coords.append(float(ycoords[0]))

            #read in z coordinate
            zrow = datarow+3
            zcoordsrow = data[zrow]
            zcoords = getcoords(zcoordsrow)
            coords.append(float(zcoords[0]))

        
        #6/9/2020: add exception for case where it splits x^-10 into x, 10
        if len(coords) > 3:
            coords_orig = coords
            coords = []
            for i in range(len(coords_orig)):
                if ( (coords_orig[i]).is_integer() == 0 or coords_orig[i] == 0): #if it's a positive integer, assume it is an exponent (this if statement is reverse of that)
                    coords.append(coords_orig[i])
                else:
                    coords[len(coords)-1] = float(coords_orig[i-1])*(10**(-int(coords_orig[i])))
            coords = np.array(coords)
        return coords

    def findWordIndices(word, data):
        wordlen = len(word)
        wordind = []
        for i in range(len(data)):
            currline = data[i]
            if (word in currline):
                wordind.append(i)
        return wordind

    #find which lines contain the relevant coordinates & check if there are >1 VOI's
    #Edit 7/20/2020: Need to include the '<' in search name because if you don't,
    #It will look for omit in the wrong place for xml files in which x, y, and z
    #coordinates each have their own line
    centerinds = findWordIndices('<Center',data)
    widthinds = findWordIndices('<WidthHalfLength',data)
    heightinds = findWordIndices('<HeightHalfLength',data)
    depthinds = findWordIndices('<DepthHalfLength',data)


    #First step: Do coordinate conversion from LPS to RCS for actual VOI
    roicenter = getXYZcoords(data,centerinds[0])
    #convert center from LPS to RAS
    roicenter[0] = -1*roicenter[0]
    roicenter[1] = -1*roicenter[1]
    

    roiradius = np.zeros((3,1))
    roiradius.astype(float)
    widthcoords = getXYZcoords(data,widthinds[0])
    roiradius[0] = widthcoords[0]

    heightcoords = getXYZcoords(data,heightinds[0])
    roiradius[1] = heightcoords[1]

    depthcoords = getXYZcoords(data,depthinds[0])
    roiradius[2] = depthcoords[2]

    #find which index of coordinates array has abs value >= 0.5, since this is what should be used for radius coordinates
    def findNonZeroInd(coords):
        for ind in range(len(coords)):
            if(round(coords[ind]) != 0):
               break
        return ind

    #Edit 7/22/2020: How to adjust radius values when x,y,z values of radius are not stored in the order width,height,depth
    #1. If widthcoords[0] or heightcoords[1] or depthcoords[2] contains a value that is approximately 0,
    #check other elements of these arrays and adjust radius array accordingly
    if( round(widthcoords[0]) == 0 or round(heightcoords[1]) == 0 or round(depthcoords[2]) == 0 ):
        #find which index of widthcoords is supposed to be used and update radius accordingly
        w_ind = findNonZeroInd(widthcoords)
        roiradius[w_ind] = widthcoords[w_ind]

        #find which index of heightcoords is supposed to be used and update radius accordingly
        h_ind = findNonZeroInd(heightcoords)
        roiradius[h_ind] = heightcoords[h_ind]

        #find which index of heightcoords is supposed to be used and update radius accordingly
        d_ind = findNonZeroInd(depthcoords)
        roiradius[d_ind] = depthcoords[d_ind]




    #Second step: If there are omit regions, get RAS coordinates for each of those
    if(len(centerinds)>1):
        ocenters = np.zeros((len(centerinds)-1,3))
        oradii = np.zeros((len(centerinds)-1,3))
          
        for i in range(1,len(centerinds)): #loop through all omit regions
            #center and half lengths in all directions for omit region of current iteration
            ocenters[i-1,:] = getXYZcoords(data,centerinds[i])
            ocenters[i-1,0] = -1*ocenters[i-1,0]
            ocenters[i-1,1] = -1*ocenters[i-1,1]

            owidthcoords = getXYZcoords(data,widthinds[i])
            oheightcoords = getXYZcoords(data,heightinds[i])
            odepthcoords = getXYZcoords(data,depthinds[i])

            oradii[i-1,0] = owidthcoords[0]
            oradii[i-1,1] = oheightcoords[1]
            oradii[i-1,2] = odepthcoords[2]
    else:
        ocenters = np.zeros((2,2))
        ocenters[0,0] = -1
        
        oradii = np.zeros((2,2))
        oradii[0,0] = -1

    return roicenter, roiradius, ocenters, oradii



def makeVOIMask(a,xs, xf, ys, yf, zs, zf, oxs, oxf, oys, oyf, ozs, ozf):
    #a is numpy array of pre-contrast image

    voi_mask = np.zeros((a.shape[0],a.shape[1],a.shape[2]))

    #give mask value 1 for all voxels in VOI
    voi_mask[xs:xf,ys:yf,zs:zf] = 1

    #check if there are any omit regions. if there aren't, oxs=oxf=oys=oyf=ozs=ozf=-1
    if (oxs[0] != -1):
        for i in range(len(oxs)):
            #give mask value 0 for all voxels in ith omit region
            voi_mask[int(oxs[i-1]):int(oxf[i-1]),int(oys[i-1]):int(oyf[i-1]),int(ozs[i-1]):int(ozf[i-1])] = 0

    return voi_mask



