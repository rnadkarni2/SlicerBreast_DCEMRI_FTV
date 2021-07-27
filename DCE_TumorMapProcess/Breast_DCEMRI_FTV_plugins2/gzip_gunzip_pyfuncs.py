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
#Script that contains all functions related to gunzip
#of compressed DICOM images.

import slicer
import os
import shutil
from distutils.dir_util import copy_tree

def extractGZ(exampath):
    
    sevenzipdir = r'C:\Program Files\7-Zip'
    if(os.path.isdir(sevenzipdir) == 0):
        slicer.util.confirmOkCancelDisplay("Error. Please decompress all files in exam directory, then try running module again.","Compressed DICOMs Error")
        return
    
    #return list of folders in exampath
    folders = [directory for directory in os.listdir(exampath) if os.path.isdir(exampath + "\\" + directory)]

    for i in range(len(folders)):
        #Edit 2/4/2021: Only bother trying to check the folder if its name is a digit,
        #because only these folders contain DICOMs.
        if(folders[i].isdigit()):
            curr_path = exampath + "\\" + folders[i]

            files = os.listdir(curr_path)
            files = [f for f in files if os.path.isfile(os.path.join(curr_path,f))] #2/8/2021: Add this step to exclude folders inside the exam folder
                                                                                    #and only take the DICOM files in that folder.
            files = sorted(files)

            #Edit 7/16/2020: Add exception for empty folders
            #Edit 11/3/2020: Add file 2 because sometimes first one is not gunzipped but rest are
            try:
                file1 = files[0]
                file2 = files[1]
                fileend = files[len(files)-1]
            except:
                file1 = ''
                file2 = ''
                fileend = ''

            #Edit 7/16/2020: create gunzipped folder inside exam folder, and in there save folders using original name
            new_dir = exampath + "\\" + "gunzipped"
            #make a gunzipped folder if necessary
            if(os.path.isdir(new_dir) == 0):
                os.mkdir(new_dir)
                
            #If folder contains gzipped DICOMs, create a new folder & add gunzipped versions of first and last image to it
            #Edit 11/3/2020: Add file2 to this if statement
            #Edit 1/26/2021: Incorporate folders in which DICOMs don't have .dcm or .DCM extension into this
            #if( ( ( file1.endswith('.dcm.gz') or file2.endswith('.dcm.gz') ) and fileend.endswith('.dcm.gz') ) or ( ( file1.endswith('.DCM.gz') or file2.endswith('.DCM.gz') ) and fileend.endswith('.DCM.gz') ) or ( ( file1.endswith('.gz') or file2.endswith('.gz') ) and fileend.endswith('.gz') )):
            #Edit 2/4/2021: new version of this if statement that also incorporates folders that have regular non-gzipped DICOMs
            if( file1.endswith('.gz') or file1.endswith('.dcm') or file1.endswith('.DCM') or file1.isdigit() ):

                new_path = new_dir + "\\" + folders[i] #save folder with original name, but inside of new gunzipped folder in exam folder

                cmd_base = r'"C:\Program Files\7-Zip\7z" x '

                #Edit 11/3/2020: if file 1 is gzipped, use the 7zip command to move it.
                #If it's just a regular DICOM, use shutil.copyfile to just copy/paste the DICOM
                if( file1.endswith('.gz') ):
                    cmd1_7z = cmd_base + curr_path + '\\' + file1 + ' -o' + new_path   #command for file 1 uncompressed
                    os.system(cmd1_7z) #execute above command
                else:
                    #2/4/2021: Edit this part to make more sense with what I have in mind now
                        #First, try seeing if file2 can be gunzipped to destination
                    if(file2.endswith('.gz')):
                        cmd2_7z = cmd_base + curr_path + '\\' + file2 + ' -o' + new_path   #command for file 1 uncompressed
                        os.system(cmd2_7z) #execute above command
                        #If file2 is also doesn't need to be gunzipped, just copy file1 to new destination without gunzipping
                    else:
                        src = curr_path + '\\' + file1
                        os.mkdir(new_path)
                        dst = new_path + '\\' + file1
                        shutil.copyfile(src,dst)

                #only run 2nd gunzip command if more than 1 file in folder
                if (len(files) > 1):
                    #Edit 2/4/2021: Also edit this part to only use gunzip if it makes sense based on file extension.
                    if(fileend.endswith('.gz')):
                        cmdend_7z = cmd_base + curr_path + '\\' + fileend + ' -o' + new_path   #command for file end uncompressed
                        os.system(cmdend_7z) #execute above command
                    else:
                        endsrc = curr_path + '\\' + fileend
                        enddst = new_path + '\\' + fileend
                        shutil.copyfile(endsrc,enddst)
    
    return



            
    


def deleteGunzipped(exampath):
    #Edit 7/16/2020: If there is a gunzipped folder inside the exampath, delete it
    gunzip_path = exampath + "\\" + "gunzipped"

    if(os.path.isdir(gunzip_path) == 1):
        shutil.rmtree(gunzip_path)

    return




#function to unzip all DCM files for relevant DCE folders, instead of just first and last DCMs
#Edit: Function only unzips all DCMs for the one DCE folder given by function inputs
def gunzipAllFilesDCE(exampath,dce_folder):

    def makeOrigAndDestPaths(exampath,dce_folder,dcegzipped_paths,dcegunzip_dest_paths):
        dce_orig_path = exampath + "\\" + dce_folder
        dcegzipped_paths.append(dce_orig_path)
        
        dce_dest_path = exampath + "\\" + "gunzipped" + "\\" + dce_folder
        dcegunzip_dest_paths.append(dce_dest_path)

        return dcegzipped_paths, dcegunzip_dest_paths
        
    dcegzipped_paths = []
    dcegunzip_dest_paths = []

    dcegzipped_paths, dcegunzip_dest_paths = makeOrigAndDestPaths(exampath,dce_folder,dcegzipped_paths,dcegunzip_dest_paths)

    #Edit 7/20/2020: First, remove existing gunzipped DICOMs from DCE folder so that user isn't prompted "replace or no"
    #Outer 'if' is extra (unnecessary) safeguard to prevent deleting original dicom or gzipped dicom images
    if('gunzipped' in dcegunzip_dest_paths[0]):
        guzfiles = os.listdir(dcegunzip_dest_paths[0])
        for file in guzfiles:
            fullfile = dcegunzip_dest_paths[0] + "\\" + file #need to use full path to file
            os.remove(fullfile)

    #Edit 2/4/2021: Edit this to accomodate cases where not all of the DICOMs in the source folder are gzipped
    files = os.listdir(dcegzipped_paths[0])
    filesdcmgz = [f for f in files if f.endswith('.gz')]
    filesdcm = [f for f in files if (f.endswith('.dcm') or f.endswith('.DCM') or f.isdigit())]

    #Using single command to unzip all DCMs in the source folder
    #Edit 2/4/2021: Only do this if all of the DICOMs are gzipped
    if(len(filesdcmgz) == len(files)):
        cmd_base = r'"C:\Program Files\7-Zip\7z" x '
        cmd_7z = cmd_base + dcegzipped_paths[0] + '\*' + ' -o' + dcegunzip_dest_paths[0]
        os.system(cmd_7z) #execute above command
    else:
        #Edit 2/4/2021: If all the DICOMs are not gzipped, use copy_tree
        if(len(filesdcm) == len(files)):
            copy_tree(dcegzipped_paths[0], dcegunzip_dest_paths[0])
        #Edit 2/4/2021: If some are gzipped and others are not, use a loop to gzip the ones
        #that are not gzipped, then move all to gunzipped destination with gunzip all in folder command.
        else:
            gzcmd_base = r'"C:\Program Files\7-Zip\7z" a -tgzip '
            #Edit 2/4/2021: if file is not gzipped, gzip it and delete the original copy
            for f_i in files:
                if('.gz' not in f_i):
                    fullfile = dcegzipped_paths[0] + '\\' + f_i 
                    gzcmd = gzcmd_base + fullfile + '.gz' + ' ' + fullfile
                    os.system(gzcmd) #create a gzipped copy of the file
                    os.remove(fullfile) #delete the original, non-gzipped copy of the file
            #Edit 2/4/2021: Once you've looped through the files to ensure that all are gzipped,
            #gunzip all of them to the new gunzipped directory.
            cmd_base = r'"C:\Program Files\7-Zip\7z" x '
            cmd_7z = cmd_base + dcegzipped_paths[0] + '\*' + ' -o' + dcegunzip_dest_paths[0]
            os.system(cmd_7z) #execute above command
                                                           
    return


