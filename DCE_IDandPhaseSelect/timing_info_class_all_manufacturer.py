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
#Script that contains functions that populate timing structures
#These timing structures are used to identify early & late
#post-contrast images

import os
import pydicom
import dicom

class GE_TimingInfo:
    def __init__(self,img1path,numpaths):
            #img1path is full path to 1st slice of the phase
            #trying to change inputs to make this compatible with both multi-folder DCE and single-folder DCE cases
        
            try:
                header = pydicom.dcmread(img1path,stop_before_pixels = True)
            except:
                header = dicom.read_file(img1path)

            try:
                self.contenttime = header[0x8,0x33].value
            except:
                self.contenttime = '000000'
                
            #content time is a good field to use to identify the phase timing for UAB as well
            #For UAB, field that tells # of slices per phase for such single folder DCE GE exams is Locations in Acquisition (21,104f)
            #Still have to check if its the same for Loyola

            try:
                self.locsinacq = header[0x21,0x104f].value
            except: #If field absent, assume each phase has its own folder and populate this field accordingly
                self.locsinacq = numpaths
            
            try:
                self.temppos = int(header[0x20,0x100].value)
                self.numtemp = int(header[0x20,0x105].value)
            except: #may have to do this for some pre-contrast cases
                self.temppos = 0
                self.numtemp = 0

            try:
                #First, try temporal resolution
                self.tempressec = header.TemporalResolution
                #convert from usec or msec to sec
                while self.tempressec > 1000:
                    self.tempressec = self.tempressec/1000
            except:
                try:
                    #Second, try acquisition duration
                    self.acqdur = header[0x19,0x105a].value #using this instead of temporal resolution header field
                    self.tempressec = self.acqdur/(self.numtemp * (10**6)) #temporal resolution in seconds is acquisition duration converted from usec to sec then divided by number of temporal positions
                except: #finally, if neither work set temporal resolution to 0 (hopefully doesn't happen for 2nd-final post-contrast
                    self.tempressec = 0
                
            try:
                self.trigtime = float(header[0x18,0x1060].value) #for ucsd image, this has same value as temporal resolution
            except:
                self.trigtime = 0

            #Edit 2/1/2021: Include Acquisition Time
            try:
                self.acqtime = header.AcquisitionTime
            except:
                self.acqtime = '000000'

            self.studydate = header.StudyDate

def getGETimingAllFolders(phaseslc1paths):
    ge_timing_all = []
    for i in range(len(phaseslc1paths)):
        curr_ge_timing = GE_TimingInfo(phaseslc1paths[i],len(phaseslc1paths))
        ge_timing_all.append(curr_ge_timing)
    return ge_timing_all




class SiemensTimingInfo:
    def __init__(self,img1path):
        try:
            header = pydicom.dcmread(img1path,stop_before_pixels = True)
        except:
            header = dicom.read_file(img1path)

        try:
            self.contenttime = header[0x8,0x33].value #Edit 1/20/2021: See if rounding this to nearest whole number prevents error for 34306 v20
        except:
            self.contenttime = '000000'

        self.studydate = header.StudyDate

        try:
            self.temppos = int(header[0x20,0x100].value)
            self.numtemp = int(header[0x20,0x105].value)
        except: #may have to do this for some pre-contrast images
            self.temppos = 0
            self.numtemp = 0
            
        try: #first, try using Temporal Resolution header field
            self.tempressec = float(header.TemporalResolution)
        except:
            try: #second, try using Scan Duration header field
                self.tempressec = float(header[0x29,0x1010].value)
            except: #third, set to 0 if neither work (hopefully only needed for some pre-contrast)
                self.tempressec = 0
            
        #convert from usec or msec to sec
        while self.tempressec >= 1000:
            self.tempressec = self.tempressec/1000

        #if temporal resolution is still > 150 after conversion to seconds, divide by number of post-contrast phases
        if self.tempressec > 150:
            if self.numtemp > 0:
                self.tempressec = self.tempressec/self.numtemp

        try:
            self.trigtime = float(header[0x18,0x1060].value) #for ucsd image, this has same value as temporal resolution
        except:
            self.trigtime = 0

        #Edit 2/1/2021: Include Acquisition Time
        try:
            self.acqtime = header.AcquisitionTime
        except:
            self.acqtime = '000000'



def getSiemensTimingAllFolders(phaseslc1paths):
    siemens_timing_all = []
    for i in range(len(phaseslc1paths)):
        curr_siemens_timing = SiemensTimingInfo(phaseslc1paths[i])
        siemens_timing_all.append(curr_siemens_timing)
    return siemens_timing_all



        
class PhilipsTimingInfo:
    def __init__(self,img1path):
        try:
            header = pydicom.dcmread(img1path,stop_before_pixels = True)
        except:
            header = dicom.read_file(img1path)

        self.studydate = header.StudyDate
        
        try:
            self.temppos = int(header[0x20,0x100].value)
            self.numtemp = int(header[0x20,0x105].value)
        except: #may have to do this for some pre-contrast images
            self.temppos = 0
            self.numtemp = 0

        try: #First, try using Temporal Resolution field for tempressec
            self.tempressec = float(header.TemporalResolution)
            #convert from usec or msec to sec
            while self.tempressec >= 1000:
                self.tempressec = self.tempressec/1000
        except: 
            try: #Second, try using Scan Duration field for tempressec
                self.tempressec = float(header[0x2005,0x1033].value)
                if self.numtemp > 0:
                    self.tempressec = self.tempressec/self.numtemp
                    #add this unit conversion part just in case it is necessary ----- convert from usec or msec to sec
                    while self.tempressec >= 1000:
                        self.tempressec = self.tempressec/1000

                else:
                    self.tempressec = 0
            except: #if tempressec is set to 0, it will later be set to difference between consecutive trigger times in chooseEarlyLate code
                self.tempressec = 0

        try:
            self.trigtime = header[0x18,0x1060].value/1000
        except:
            self.trigtime = 0

        #Edit 6/14/21: Need to add content time and acquisition time to this code.
        try:
            self.contenttime = header[0x8,0x33].value #Edit 1/20/2021: See if rounding this to nearest whole number prevents error for 34306 v20
        except:
            self.contenttime = '000000'

        try:
            self.acqtime = header.AcquisitionTime
        except:
            self.acqtime = '000000'


def getPhilipsTimingAllFolders(phaseslc1paths):
    philips_timing_all = []
    for i in range(len(phaseslc1paths)):
        curr_philips_timing = PhilipsTimingInfo(phaseslc1paths[i])
        philips_timing_all.append(curr_philips_timing)
    return philips_timing_all
        

        
            

