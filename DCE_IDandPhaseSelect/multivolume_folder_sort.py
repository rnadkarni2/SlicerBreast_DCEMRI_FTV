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
#Script that contains functions for sorting DICOM slices in a
#multivolume DCE folder in the correct order

#Edit 7/24/2020: Using content time to sort into different phases, instead of temporal position

import os
import pydicom
import numpy as np


#Class of image info object. Each object stores the name of the DICOM file
#as well as some header info

#For Philips, you use Temporal Position Identifier and Number of Temporal Positions
#For GE, you use Content Time for 1st sort instead
class DicomFileInfo:
    def __init__(self,dcepath,file):
        self.filename = file
        
        imgpath = dcepath + "\\" + file
        hdr = pydicom.dcmread(imgpath,stop_before_pixels = True)
        try:
            self.numtemp = int(hdr[0x20,0x105].value)
        except:
            self.numtemp = 0
        try:
            self.temppos = int(hdr[0x20,0x100].value)
        except:
            #if field is absent, assume pre-contrast if Philips exam
            self.temppos = 0

        #6/13/2021: Try to implement David's backup option for when
        #Slice Location header field is missing.
        try:
            self.slcloc = float(hdr[0x20,0x1041].value)
        except:
            img_orient = hdr[0x20,0x37].value
            plane_cos = np.cross(img_orient[0:3],img_orient[3:])
            img_pos = hdr[0x20,0x32].value
            self.slcloc = np.dot(plane_cos,img_pos)

        try:
            self.ctime = float(hdr.ContentTime)
        except:
            self.ctime = 0
        
        try:
            self.trigtime = float(hdr.TriggerTime)
        except:
            self.trigtime = 0

        #Edit 2/1/2021: Acquisition Time is the 3rd option when Content Time and Trigger Time don't work
        try:
            self.acqtime = float(hdr.AcquisitionTime)
        except:
            self.acqtime = 0

        #Edit 2/26/2021: Add Images in Acquisition
        try:
            self.im_in_acq = int(hdr[0x20,0x1002].value)
        except:
            self.im_in_acq = 0

        #Edit 4/2/2021: Add manufacturer, because trigger time should be default
        #phase timing field for Philips, and content time should be default for
        #other manufacturers.
        self.manufacturer = hdr.Manufacturer

        

        


#dcepath = r"\\researchfiles.radiology.ucsf.edu\birp_perm2\ispy_2019\4460_OHSC\13432\20190910_v10\E20190910\701"
#it's faster to use os.listdir that sitk method of retrieving filenames

def sortDicoms(dcepath):
    print("Starting multivolume DICOM sort........................")
    dicom_names = os.listdir(dcepath)
    dicom_names = sorted(dicom_names)

    totslices = len(dicom_names)

    #Use dicom names to create unsorted list of DicomFileInfo class objects
    all_dcms_info = []
    contenttimesall = []
    trigtimesall = []
    acqtimesall = []
    for dicom_name in dicom_names:
        curr_dcm_info = DicomFileInfo(dcepath,dicom_name)
        all_dcms_info.append(curr_dcm_info)
        contenttimesall.append(curr_dcm_info.ctime)
        trigtimesall.append(curr_dcm_info.trigtime)
        acqtimesall.append(curr_dcm_info.acqtime)

    im_in_acq = curr_dcm_info.im_in_acq
    #all_dcms_info = sorted(all_dcms_info, key = lambda s: s.temppos) #sort by phase

    #Take DICOM file info structures for first and last slices. If both have numtemp greater than 0, then nslice = len(dicom_names)/numtemp
    #If one of them has numtemp == 0, real numtemp = numtemp + 1 and nslice is calculated with this numtemp
##    dcm1info = all_dcms_info[0]
##    dcmendinfo = all_dcms_info[len(all_dcms_info)-1]
##
##    if(dcm1info.numtemp > 0 and dcmendinfo.numtemp > 0):
##        numtemp = int(dcmendinfo.numtemp)
##    if(dcm1info.numtemp == 0 or dcmendinfo.numtemp == 0):
##        numtemp = int(dcmendinfo.numtemp) + 1

    #Edit 7/24/2020: use number of unique content time values as numtemp
    contenttimesall = np.array(contenttimesall)

    print("content times all")
    print(contenttimesall)
        

    print("content times all after removing small differences between adjacent elements")
    print(contenttimesall)
    
    contenttimesunique = np.unique(contenttimesall)
    print("content times unique")
    print(contenttimesunique)
    #Edit 1/28/2021: Use trigger time when < 3 unique content times, not just when there is only 1 content time.
    #Edit 2/1/2021: Use very small range of content times (< 1 min) as a condition for using other fields for timing
    #Edit 2/2/2021: Add these binary variables to synchronize choice of timing variable between multivolume_folder_sort and choose_early_late
    ctime_forphases = 0
    ttime_forphases = 0
    atime_forphases = 0

    try:
        nphase_forcheck = int(totslices/im_in_acq)
    except:
        nphase_forcheck = 0 #have to do this if unable to retrieve im_in_acq from header
        
    maxdiff = (np.amax(contenttimesunique)-np.amin(contenttimesunique))
    print("number of phases")
    print(nphase_forcheck)
    print("content times")
    print(contenttimesunique)
    print("max diff in content times")
    print(maxdiff)
    
    if('Philips' in curr_dcm_info.manufacturer or 'PHILIPS' in curr_dcm_info.manufacturer):
        ttime_forphases = 1 #4/2/2021: For now, assume Philips always uses trigger time
        trigtimesall = np.array(trigtimesall)
        contenttimesunique = np.unique(trigtimesall)
        #Edit 6/14/2021: Even for Philips, need to check if trigger time actually
        #works and use other timing variable if it doesn't
        if(len(contenttimesunique) <= 2):
            contenttimesunique = np.unique(contenttimesall)
            
            ctime_forphases = 1
            ttime_forphases = 0
            if(len(contenttimesunique) <= 2):
                print("Using acquisition time to sort multivolume folder slices by DCE phase")
                contenttimesunique = np.unique(acqtimesall)
                atime_forphases = 1
                ctime_forphases = 0
            else:
                print("Using content time to sort multivolume folder slices by DCE phase")

        else:
            print("Using trigger time to sort multivolume folder slices by DCE phase")


        
    #Edit 4/2/2021: Only default to content time for phase sorting for Siemens and GE
    else:
        #Edit 2/26/2021: Added new check --> use totslices and im_in_acq to make sure that number of content times matches number of phases.
        #Edit 3/23/2021: Only use content times if totslices is divisible by # of content times
        if( len(contenttimesunique) < 3 or maxdiff<59 or (im_in_acq>0 and nphase_forcheck > 1 and len(contenttimesunique) < nphase_forcheck ) or (totslices%len(contenttimesunique) != 0) ):
            trigtimesall = np.array(trigtimesall)
            #Edit 2/1/2021: If trigger times field is empty, use Acquisition Times
            if(np.amax(trigtimesall) > 0):
                contenttimesunique = np.unique(trigtimesall)
                print("Using trigger time to sort multivolume folder slices by DCE phase")
                ttime_forphases = 1
                #all_dcms_info = sorted(all_dcms_info, key = lambda s: s.trigtime) #Edit 7/24/2020: use trigger times to sort by phase
            else:
                contenttimesunique = np.unique(acqtimesall)
                print("Using acquisition time to sort multivolume folder slices by DCE phase")
                atime_forphases = 1
                #all_dcms_info = sorted(all_dcms_info, key = lambda s: s.acqtime) #Edit 2/1/2021: use trigger times to sort by phase

        else:
            print("Using content time to sort multivolume folder slices by DCE phase")
            ctime_forphases = 1
            #all_dcms_info = sorted(all_dcms_info, key = lambda s: s.ctime) #Edit 7/24/2020: use content times to sort by phase

    #3/22/2021: Try to get rid of small timing differences (1 second) before sorting slices
    dlt_times1 = []
    for t in range(1,len(contenttimesunique)):
        currtime = contenttimesunique[t]
        #need 2nd case in this if statement because haven't converted from hhmmss to sec
        if(currtime - contenttimesunique[t-1] < 2 or (currtime - contenttimesunique[t-1] == 41 and str(int(contenttimesunique[t-1])).endswith('59')) ):
            contenttimesunique[t] = contenttimesunique[t-1]
            dlt_times1.append(currtime)
    contenttimesunique = np.unique(contenttimesunique)

    if(len(dlt_times1) > 0):
        for s in range(len(all_dcms_info)):
            if(ttime_forphases == 1):
                if(all_dcms_info[s].trigtime in dlt_times1):
                    if(str(int(all_dcms_info[s].trigtime)).endswith('00')):
                        truetime = all_dcms_info[s].trigtime - 41
                    else:
                        truetime = all_dcms_info[s].trigtime - 1
                    all_dcms_info[s].trigtime = truetime

            if(atime_forphases == 1):
                if(all_dcms_info[s].acqtime in dlt_times1):
                    if(str(int(all_dcms_info[s].acqtime)).endswith('00')):
                        truetime = all_dcms_info[s].acqtime - 41
                    else:
                        truetime = all_dcms_info[s].acqtime - 1
                    all_dcms_info[s].acqtime = truetime
                    
            if(ctime_forphases == 1):
                if(all_dcms_info[s].ctime in dlt_times1):
                    if(str(int(all_dcms_info[s].ctime)).endswith('00')):
                        truetime = all_dcms_info[s].ctime - 41
                    else:
                        truetime = all_dcms_info[s].ctime - 1
                    all_dcms_info[s].ctime = truetime
        
    #3/22/21: if timing array is longer than nphase, use the first nphase elements only
    if(nphase_forcheck > 1 and len(contenttimesunique) > nphase_forcheck):
        for v in range(len(all_dcms_info)):
            if(ttime_forphases == 1):
                if(all_dcms_info[v].trigtime > contenttimesunique[nphase_forcheck-1]):
                    all_dcms_info[v].trigtime = contenttimesunique[nphase_forcheck-1]
            if(atime_forphases == 1):
                if(all_dcms_info[v].acqtime > contenttimesunique[nphase_forcheck-1]):
                    all_dcms_info[v].acqtime = contenttimesunique[nphase_forcheck-1]                
            if(ctime_forphases == 1):
                if(all_dcms_info[v].ctime > contenttimesunique[nphase_forcheck-1]):
                    all_dcms_info[v].ctime = contenttimesunique[nphase_forcheck-1]                
        contenttimesunique = contenttimesunique[0:nphase_forcheck]

    #Edit 3/22/2021: only sort slices after fixes above
    if(ttime_forphases == 1):
        all_dcms_info = sorted(all_dcms_info, key = lambda s: s.trigtime) #Edit 7/24/2020: use trigger times to sort by phase
    if(atime_forphases == 1):
        all_dcms_info = sorted(all_dcms_info, key = lambda s: s.acqtime)
    if(ctime_forphases == 1):
        all_dcms_info = sorted(all_dcms_info, key = lambda s: s.ctime)
        
    numtemp = len(contenttimesunique)
    print("number of temporal positions based on timing variable")
    print(numtemp)
    nslice = int(len(dicom_names)/(numtemp)) #number of slices per phase
    print("number of slices per phase using numtemp described above")
    print(nslice)
    print("content times unique")
    print(contenttimesunique)

    #create 2D array that will store filenames in order
    w,h = nslice,numtemp #row index is phase number, column index is slice number
    fsort = [['' for x in range(w)] for y in range(h)]
    #fill this array row by row
    for r in range(numtemp):
        phase_dcms_info = all_dcms_info[r*nslice:(r+1)*nslice] #all dicom objects for current phase
        phase_dcms_info = sorted(phase_dcms_info, key = lambda s: s.slcloc, reverse = True) #sort all slices for current phase in reverse order of slice location

        for c in range(nslice):
            slc_info = phase_dcms_info[c] #for current phase, take header object for cth slice
            fsort[r][c] = slc_info.filename #save filename from header object into appropriate position in string array according to phase and slice order

    #This works! Although using slice location appears to sort images in reverse order-- use reverse = True for this 1st sort
    #Apparently, this makes it correct order for UChic exams and reverse order for OHSC exams -- don't bother doing anything about this for now
    #print(fsort[0][:])
    print(".......................... Finished multivolume DICOM sort")
    return fsort, numtemp, nslice, ctime_forphases, ttime_forphases, atime_forphases
    


    

    



    
