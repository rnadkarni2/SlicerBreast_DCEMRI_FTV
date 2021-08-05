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
#Script that contains functions for identifying
#early and late post-contrast images



import time
import pydicom
import dicom
import os
import numpy as np

import Breast_DCEMRI_FTV_plugins1
from Breast_DCEMRI_FTV_plugins1 import timing_info_class_all_manufacturer #function for populating structures of timing info from DICOM headers
from Breast_DCEMRI_FTV_plugins1 import multivolume_folder_sort #function for sorting slices in correct order for a multivolume folder
from Breast_DCEMRI_FTV_plugins1 import gzip_gunzip_pyfuncs

#Function that tells you which post-contrast image # is early or late
def findClosestTime(imgtimes,time):
    timediffs = np.abs(imgtimes-time)
    ind = (np.where(timediffs==np.min(timediffs)))
    print("ind is ")
    print(ind)

    #try-except statements to prevent error at this step
    try:
        postContrastNum = 1 + int(ind[0])
    except:
        postContrastNum = 1 + ind[0][0]

    return postContrastNum


#Function that tells you how many mm:ss after content time of 1st post-contrast early or late post-contrast is
#Output of this function is only used in report PDF
def getTimeMMSS(timediff):
    timediffmm = int(timediff/60)
    timediffss = timediff%60
    try:
        timediffss = int(round(timediffss))
    except:
        timediffss = int(round(timediffss[0]))

    return timediffmm, timediffss


#Function for choosing early and late post-contrast phases for
#GE and Siemens exams
#7/6/2021: make number of seconds post-contrast
#for early and late function inputs, instead of fixed values of 150 and 450
def chooseEarlyLateGE_Siemens(exampath,dce_folders,manufacturer,earlyadd, lateadd):
    #for gunzip all DCMs, need to use original exampath as input
    orig_exampath = exampath[:-10] #remove \gunzipped from the end to get original exampath

    #Edit 7/17/2020: If 'gunzipped' in exampath, gunzip all DCMs for precontrast folder (or in some cases, DCE multivolume folder)
    if('gunzipped' in exampath):
        gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[0]))


    #use one slice from each DCE folder
    if (len(dce_folders) > 2):
        nslice = 0 #setting this to 0 for multi-folder DCE makes this compatible with other functions used for Slicer
        fsort = 0 #setting this to 0 for multi-folder DCE makes this compatible with other functions used for Slicer
        phaseslc1paths = []
        for fnum in range(len(dce_folders)):
            curr_path = os.path.join(exampath,str(dce_folders[fnum]))
            files = [f for f in os.listdir(curr_path) if f.endswith('.dcm')]
            FILES = [f for f in os.listdir(curr_path) if f.endswith('.DCM')]
            files_noext = [f for f in os.listdir(curr_path) if f.isdigit()] #edit 1/26/2021: In some folders, there is a series of DICOM images
                                                                            #with no .dcm or .DCM extension, but all the files have a number as the filename.

            if len(files) > 0:
                files = sorted(files)
                curr_img_path = os.path.join(curr_path,files[0])
            if len(FILES) > 0:
                FILES = sorted(FILES)
                curr_img_path = os.path.join(curr_path,FILES[0])
            if len(files_noext)>0:
                curr_img_path = os.path.join(curr_path,files_noext[0])


            phaseslc1paths.append(curr_img_path)

    #use one slice from each phase in same folder if there is only one DCE folder
    if (len(dce_folders) == 1):
        dcepath = os.path.join(exampath,str(dce_folders[0]))
        fsort, numtemp, nslice, ctime_forphases, ttime_forphases, atime_forphases = multivolume_folder_sort.sortDicoms(dcepath)
        phaseslc1paths = []
        for p in range(len(fsort)):
            curr_img_path = os.path.join(dcepath,fsort[p][0])
            phaseslc1paths.append(curr_img_path)
        print("Phases slice 1 paths")
        print(phaseslc1paths)

    #Edit 7/17/2020: use one slice from precontrast folder and one slice from each phase in multivolume postcontrast folder
    if(len(dce_folders) == 2):
        #Edit 7/17/2020: If 'gunzipped' in exampath, gunzip all DCMs for multivolume postcontrast folder
        if('gunzipped' in exampath):
            gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[1]))

        postpath = os.path.join(exampath,str(dce_folders[1]))
        fsort, numtemp, nslice, ctime_forphases, ttime_forphases, atime_forphases = multivolume_folder_sort.sortDicoms(postpath)

        prepath = os.path.join(exampath,str(dce_folders[0]))
        prefiles = os.listdir(prepath)
        prefiles = sorted(prefiles)
        preslc1path = os.path.join(prepath,prefiles[0])

        phaseslc1paths = []
        phaseslc1paths.append(preslc1path)
        for p in range(len(fsort)):
            curr_img_path = os.path.join(postpath,fsort[p][0])
            phaseslc1paths.append(curr_img_path)

    #Create array of Timing Structures
    if ('GE' in manufacturer):
        timing_all = timing_info_class_all_manufacturer.getGETimingAllFolders(phaseslc1paths)
    if ('Siemens' in manufacturer or 'SIEMENS' in manufacturer):
        timing_all = timing_info_class_all_manufacturer.getSiemensTimingAllFolders(phaseslc1paths)

    #initialize array of times (in seconds) for all post-contrast images
    contenttimes = np.zeros((len(timing_all)-1,1))
    trigtimes = np.zeros((len(timing_all)-1,1))
    acqtimes = np.zeros((len(timing_all)-1,1)) #Edit 2/1/2021: Acquisition Time
    use_trig_time = 0
    use_acq_time = 0
    print("length timing all")
    print(len(timing_all))

    #6/24/21: If you change timing variable from content time to other and
    #tempres is currently contenttime(2nd post) - contenttime(1st post), you
    #have to change tempres to match the new timing variable
    tempres_subtract = 0

    #skip pre-contrast image, and return array of post-contrast image content times, in seconds
    for i in range(1,len(timing_all)):
        curr_timing = timing_all[i]

        ctime = curr_timing.contenttime
        #Edit 6/27/2021: If ctime has less than 6 digits, it is not a valid
        #hh:mm:ss time, so then make ctimesec 0 so that one of the other
        #phase timing variables will be used.
        print("ctime")
        print(ctime)
        print(ctime[0])
        #Edit 7/2/2021: Getting rid of this 6 digits check because it
        #causes problems for exams where hours value is < 10.
        try:
            ctimesec = 3600*int(ctime[0:2]) + 60*int(ctime[2:4]) + float(ctime[4:])
        except: #set ctimesec to 0 if above gives error due to ctime < 5 digits
            ctimesec = 0
        contenttimes[i-1,0] = ctimesec
        trigtimes[i-1,0] = curr_timing.trigtime/1000 #For UKCC, must use this in place of content time

        acqtime = curr_timing.acqtime
        acqtimesec = 3600*int(acqtime[0:2]) + 60*int(acqtime[2:4]) + float(acqtime[4:]) #For UMinn 23929 v10, must use this in place of content time.
        acqtimes[i-1,0] = acqtimesec

        #Get temporal resolution (in seconds) from 2nd post-contrast image. Also read Study Date from this image.
        if i == 2:
            tempres = curr_timing.tempressec
            studydate = curr_timing.studydate

            #Edit 6/12/2020: Like Philips, make tempres equal to difference between trigger (content) times if necessary
            #Edit 3/22/2021: Also use difference between 2nd and 1st post-contrast
            #timing as tempres if tempres from header is too large (> 120 sec).
            if ( tempres < 30 or tempres > 120): #Edit: use small fixed value threshold instead of actually setting it to difference between trigger/content times
                tempres = contenttimes[i-1,0] - contenttimes[i-2,0]
                #For UKCC, use trigger time because all content times are the same
                if(tempres == 0):
                    use_trig_time = 1
                    tempres = trigtimes[i-1,0] - trigtimes[i-2,0]
                tempres_subtract = 1

    print("Content times just before trigger time check")
    print(contenttimes)

    #Edit 2/2/2021: If multivolume_folder_sort was used, make sure choice of timing variable is same as the one chosen in that function.
    if(len(dce_folders) == 1 or len(dce_folders) == 2):
        use_trig_time = ttime_forphases
        use_acq_time = atime_forphases
    else:
        #For some exams, like UKCC, must use trigger times in place of content times
        #Edit 2/1/2021: Add scenarios for using acqtimes

        #Edit 2/8/2021: Added separation of 59 seconds between min and max content times as another condition
        #for choosing timing variable.
        try:
            if(contenttimes[1,0] == contenttimes[2,0] or (np.amax(contenttimes)-np.amin(contenttimes) < 59) ):
                if(np.amax(trigtimes) > 0):
                    use_trig_time = 1
                else:
                    use_acq_time = 1
        except:
            #Edit 1/21/2021: For cases like UAB exam 19701 v30 where there are only 2 post-contrast phases
            if(contenttimes[0,0] == contenttimes[1,0] or(np.amax(contenttimes)-np.amin(contenttimes) < 59) ):
                if(np.amax(trigtimes) > 0):
                    use_trig_time = 1
                else:
                    use_acq_time = 1

    #For some exams, like UKCC, must use trigger times in place of content times
    if(use_trig_time == 1):
        print("using trigger time in place of content time")
        contenttimes = trigtimes
        #6/24/21: Change tempres to match new timing variable if necessary
        if(tempres_subtract == 1):
            tempres = trigtimes[1,0] - trigtimes[0,0]

    #Edit 2/1/2021: For some exams, like UMinn 23929 v10, must use acquisition times in place of content times
    if(use_acq_time == 1):
        print("using acquisition time in place of content time")
        contenttimes = acqtimes
        #6/24/21: Change tempres to match new timing variable if necessary
        if(tempres_subtract == 1):
            tempres = acqtimes[1,0] - acqtimes[0,0]

    #3/22/2021: Try to get rid of small timing differences (1 second)
    #before doing early/late timing id
    smalldiff = 0
    for t in range(1,len(contenttimes)):
        currtime = contenttimes[t]
        #don't need 2nd case in 'if' here because times have
        #already been converted from hhmmss to sec
        if(currtime - contenttimes[t-1] < 2):
            smalldiff = 1
            contenttimes[t] = contenttimes[t-1]
    contenttimes = np.unique(contenttimes)
    if(smalldiff == 1):
        tempres = contenttimes[2]-contenttimes[1]

    print("content times")
    print(contenttimes)
    print(" ")
    print("temporal resolution")
    print(tempres)
    print(" ")

    imgtimes = contenttimes+0.5*tempres #timings used for early/late labeling

    if(len(contenttimes) > 2):
        #7/6/2021: change to user input values instead of always
        #2.5*60 and 7.5*60
        earlytime = contenttimes[0]+earlyadd #2.5*60 #early post-contrast image is closest to 2.5 min after content time of 1st post-contrast
        latetime = contenttimes[0]+lateadd #7.5*60 #early post-contrast image is closest to 2.5 min after content time of 1st post-contrast

        print("image times")
        print(imgtimes)
        print(" ")
        print("early time")
        print(earlytime)
        print(" ")
        print("late time")
        print(latetime)
        print(" ")

        earlyPostContrastNum = findClosestTime(imgtimes,earlytime)
        latePostContrastNum = findClosestTime(imgtimes,latetime)
    else:
        #Edit 1/21/2021: For cases like UAB exam 19701 v30 where there are only 2 post-contrast phases
        earlyPostContrastNum = 1
        latePostContrastNum = 2

    #Edit 2/26/2021: If obvious error such as same phase labeled early and late,
    #Redo timing ident using only numtemp and tempres
    #Edit 3/16/2021: Adding another case for use of numtemp and tempres, because if 1st post-contrast
    #is chosen as early that usually means timing identification is wrong.
    if(earlyPostContrastNum == latePostContrastNum or (earlyPostContrastNum == 1 and len(contenttimes) > 2) ):
        tempposits = np.array(range(curr_timing.numtemp))
        imgtimes_new = tempres*tempposits + 0.5*tempres
        earlytime= 2.5*60
        latetime = 7.5*60
        #To prevent error, check length of imgtimes before trying
        #to replace existing values of early and late nums
        if(len(imgtimes_new) > 2):
            earlyPostContrastNum = findClosestTime(imgtimes_new,earlytime)
            latePostContrastNum = findClosestTime(imgtimes_new,latetime)

    #If more than 1 DCE folder, gunzip all DCMs for early and late post-contrast folders
    if(len(dce_folders) > 2):
        #Edit 7/17/2020: If 'gunzipped' in exampath, gunzip all DCMs for early post-contrast folder
        if('gunzipped' in exampath):
            gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[earlyPostContrastNum]))

        #Edit 7/17/2020: If 'gunzipped' in exampath, gunzip all DCMs for late post-contrast folder
        if('gunzipped' in exampath):
            gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[latePostContrastNum]))

    #find out how many min and sec after start of 1st post-contrast early post-contrast occurs, for report
    earlydiff = imgtimes[earlyPostContrastNum-1] - contenttimes[0]
    earlydiffmm, earlydiffss = getTimeMMSS(earlydiff)

    #find out how many min and sec after start of 1st post-contrast late post-contrast occurs, for report
    latediff = imgtimes[latePostContrastNum-1] - contenttimes[0]
    latediffmm, latediffss = getTimeMMSS(latediff)

    return tempres, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss


#Function for choosing early and late post-contrast phases for
#Philips exams
#7/6/2021: make number of seconds post-contrast
#for early and late function inputs, instead of fixed values of 150 and 450
def chooseEarlyLatePhilips(exampath,dce_folders,earlyadd,lateadd):
    #for gunzip all DCMs, need to use original exampath as input
    orig_exampath = exampath[:-10] #remove \gunzipped from the end to get original exampath

    print("dce folders at start of Philips choose early/late")
    print(dce_folders)

    #Edit 7/17/2020: If 'gunzipped' in exampath, gunzip all DCMs for DCE multivolume folder
    if('gunzipped' in exampath):
        gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[0]))
        #Edit 3/29/2021: If 2 DCE folders, must also gunzip the 2nd one
        if(len(dce_folders) == 2):
            gzip_gunzip_pyfuncs.gunzipAllFilesDCE(orig_exampath,str(dce_folders[1]))

    #Using test_multivolume_folder_sort to automatically sort dicom filenames
    if(len(dce_folders) == 1):
        dcepath = os.path.join(exampath,str(int(dce_folders[0])))

    if(len(dce_folders) == 2):
        dcepath = os.path.join(exampath,str(int(dce_folders[1])))

    print("dce path")
    print(dcepath)

    fsort,numtemp,nslice, ctime_forphases, ttime_forphases, atime_forphases = multivolume_folder_sort.sortDicoms(dcepath)

    print("ttime for phases")
    print(ttime_forphases)
    print("ctime for phases")
    print(ctime_forphases)
    print("atime for phases")
    print(atime_forphases)

    #For Philips, can assume always single-folder DCE
    #3/29/2021: Brtool results for OHSC 06294 v10 showed that
    #Philips can have 2-folder DCE too. I'm editing the code
    #to allow processing for these cases.
    phaseslc1paths = []
    if(len(dce_folders) == 2):
        prepath = os.path.join(exampath,str(dce_folders[0]))
        prefiles = os.listdir(prepath)
        prefiles = sorted(prefiles)
        preslc1path = os.path.join(prepath,prefiles[0])
        phaseslc1paths.append(preslc1path)

    for i in range(len(fsort)):
        curr_img_path = os.path.join(dcepath,fsort[i][0])
        phaseslc1paths.append(curr_img_path)

    print("1st slice of each phase")
    print(phaseslc1paths)

    timing_all = timing_info_class_all_manufacturer.getPhilipsTimingAllFolders(phaseslc1paths)

    trig_times = []
    #loop to obtain unique trigger time values & store them in trig_times
    #skip precontrast, and read 1st slice for each post-contrast phase
    for j in range(1,len(timing_all)):
        curr_timing = timing_all[j]
        #Edit 6/14/2021: Need to edit this code to account for
        #fact that trigger time is not the timing variable for
        #Philips exams on rare occasions.
        if(ttime_forphases == 1):
            trig_times.append(curr_timing.trigtime)
        if(ctime_forphases == 1):
            ctime = curr_timing.contenttime
            ctime_sec = 3600*int(ctime[0:2]) + 60*int(ctime[2:4]) + float(ctime[4:])
            trig_times.append(ctime_sec)
        if(atime_forphases == 1):
            atime = curr_timing.acqtime
            atime_sec = 3600*int(ctime[0:2]) + 60*int(ctime[2:4]) + float(ctime[4:])
            trig_times.append(atime_sec)

        if j == 2:
            studydate = curr_timing.studydate
            tempres = curr_timing.tempressec

            #set tempres to difference between consecutive trig times if necessary
            if ( tempres < 30 ): #Edit: use small fixed value threshold instead of actually setting it to difference between trigger/content times
                tempres = trig_times[j-1] - trig_times[j-2]

    #sort slices in correct order based on trigger time
    trig_times = np.array(trig_times)
    print("trigger times")
    print(trig_times)

    imgtimes = trig_times + 0.5*tempres #Philips version of content time + 0.5*temporal resolution

    #7/6/2021: change to user input values instead of always
    #2.5*60 and 7.5*60
    earlytime = trig_times[0]+earlyadd #2.5*60 #early post-contrast image is closest to 2.5 min after content time of 1st post-contrast
    latetime = trig_times[0]+lateadd #7.5*60 #early post-contrast image is closest to 2.5 min after content time of 1st post-contrast

    print(imgtimes)
    print(earlytime)
    print(latetime)

    earlyPostContrastNum = findClosestTime(imgtimes,earlytime)
    print(earlyPostContrastNum)

    latePostContrastNum = findClosestTime(imgtimes,latetime)
    print(latePostContrastNum)

    #find out how many min and sec after start of 1st post-contrast early post-contrast occurs, for report
    earlydiff = imgtimes[earlyPostContrastNum-1] - trig_times[0]
    earlydiffmm, earlydiffss = getTimeMMSS(earlydiff)

    #find out how many min and sec after start of 1st post-contrast late post-contrast occurs, for report
    latediff = imgtimes[latePostContrastNum-1] - trig_times[0]
    latediffmm, latediffss = getTimeMMSS(latediff)

    return tempres, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss

