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

#This is the 1st of 2 modules in the extension FTV_process_complete.
#The purpose of this module is to load the pre-contrast, early post-contrast,
#and late post-contrast phases of the MR exam selected by the user.

#The user can either allow the module to automatically identify
#the DCE series or manually select the series that they want to
#be identifed as DCE. In addition, the user has the option to
#change target times the module uses to select early and late
#post-contrast phases.


import os
import sys
import unittest
import vtk, qt, ctk, slicer
from qt import *
from slicer.ScriptedLoadableModule import *
import logging

try:
  import numpy as np
except:
  slicer.util.pip_install('numpy')
  import numpy as np

try:
  import nibabel as nib
except:
  slicer.util.pip_install('nibabel')
  import nibabel as nib

from vtk.util.numpy_support import vtk_to_numpy
from vtk.util import numpy_support
#import scipy
#from scipy import signal
import qt
import re

try:
  import pydicom
except:
  slicer.util.pip_install('pydicom')
  import pydicom

try:
  import dicom
except:
  slicer.util.pip_install('dicom')
  import dicom

try:  
  import scipy
except:
  slicer.util.pip_install('scipy')
  import scipy
  
from scipy import signal

import math
  
import time
import datetime
import pickle
import pathlib

#imports of my functions
#7/27/2021: Start with Plugin folder import
import Breast_DCEMRI_FTV_plugins1
from Breast_DCEMRI_FTV_plugins1 import Get_header_info_all_manufacturer
from Breast_DCEMRI_FTV_plugins1 import read_DCE_images_to_numpy
from Breast_DCEMRI_FTV_plugins1 import Exam_Ident_and_timing
from Breast_DCEMRI_FTV_plugins1 import ident_gzipped_exam
from Breast_DCEMRI_FTV_plugins1 import gzip_gunzip_pyfuncs

#
# DCE_IDandPhaseSelect
#





#Code for finding pre, early, and late images and adding them to nodes
#Separate function for this now (7/28/2020) because it may have to be
#called for multiple visits

#Edit 12/7/2020: Added binary variable orig
#For original visit selected by user, this will have value of 1 indicating
#that workspace of Exam_ident_and_timing outputs should be saved.
#If user loads other visits, this will have value of 0 for those visits so that
#you don't save workspace for those visits.

def loadPreEarlyLate(exampath,visitnum,orig,dce_folders_manual,dce_ind_manual,earlyadd,lateadd):

  #6/28/2021: Make this compatible with exams that don't use
  #\\researchfiles directory structure
  if('ispy' in exampath or 'acrin' in exampath):
    #Do this so you can have MR2.5 but MR1 will not be written as MR1.0
    if(visitnum%10 == 0):
      mrnum = int(visitnum/10)
    else:
      mrnum = float(visitnum/10)
    nodevisstr = 'MR'+str(mrnum)+' ' #v10 = MR1, v20 = MR2, etc
    prenodestr = nodevisstr + 'pre-contrast'
    earlynodestr = nodevisstr + 'early post-contrast'
    latenodestr = nodevisstr + 'late post-contrast'

  else:
    nodevisstr = visitnum
    prenodestr = nodevisstr + ' pre-contrast'
    earlynodestr = nodevisstr + ' early post-contrast'
    latenodestr = nodevisstr + ' late post-contrast'

  

  #Edit 5/8/2020: Code for identifying manufacturer and then calling manufacturer-specific
  #function for identifying pre-contrast and post-contrast DCE folder(s)
    
  #Edit 6/11/2020: Call function that does DCE folder and early/late timing identification given exampath
  tempres, all_folders_info, dce_folders, dce_ind, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss = Exam_Ident_and_timing.runExamIdentAndTiming(exampath,dce_folders_manual,dce_ind_manual,earlyadd,lateadd)

  #edit 5/19/2020: stop saving .nii files, and only load the images to Slicer nodes
  #Edit 8/7/2020: Change progress bar title to include visit information
  #Edit 7/2/2021: Only create progressBar if there are DCE folders found.
  if(len(dce_folders) > 0):
    progbarttl = "Visit: " + nodevisstr
    progressBar = slicer.util.createProgressDialog(windowTitle = progbarttl)


  #Save workspace for original exam selected by user
  if(orig == 1):
    #7/9/2021: Make this path part independent of computer & directory
    #code is stored in.
    extension_path = os.path.join( os.path.dirname( __file__ ), '..' )
    wkspc_savepath = os.path.join(extension_path,'current_exam_workspace.pickle')
    with open(wkspc_savepath,'wb') as f:
      pickle.dump([tempres, all_folders_info, dce_folders, dce_ind, fsort, studydate, nslice, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss],f)

  #Edit 11/25/2020: save all outputs from exam_ident_and_timing to parameter nodes
##  tempres_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","tempres node")
##  tempres_node.SetParameter("tempres",str(tempres))
##
##  allfolders_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","all folders node")
##  allfolders_node.SetParameter("all_folders_info",all_folders_info)
##
##  dcefolders_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","dce folders node")
##  dcefolders_node.SetParameter("dce_folders",dce_folders)
##
##  dceind_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","dce indices node")
##  dceind_node.SetParameter("dce_ind",dce_ind)
##
##  fsort_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","fsort node")
##  fsort_node.SetParameter("fsort",fsort)
##
##  studydate_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","studydate node")
##  studydate_node.SetParameter("studydate",studydate)
##
##  nslice_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","nslice node")
##  nslice_node.SetParameter("nslice",nslice)
##
##  earlyPostContrastNum_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","earlyPostContrastNum node")
##  earlyPostContrastNum_node.SetParameter("earlyPostContrastNum",earlyPostContrastNum)
##
##  earlydiffmm_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","earlydiffmm node")
##  earlydiffmm_node.SetParameter("earlydiffmm",earlydiffmm)
##
##  earlydiffss_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","earlydiffss node")
##  earlydiffss_node.SetParameter("earlydiffss",earlydiffss)
##
##  latePostContrastNum_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","latePostContrastNum node")
##  latePostContrastNum_node.SetParameter("latePostContrastNum",latePostContrastNum)
##
##  latediffmm_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","latediffmm node")
##  latediffmm_node.SetParameter("latediffmm",latediffmm)
##
##  latediffss_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","latediffss node")
##  latediffss_node.SetParameter("latediffss",latediffss)
  #End edit 11/25/2020

  lbltxt = "DCE images are in folders " + str(dce_folders) + ". Early and late post-contrast phases are " + str(earlyPostContrastNum) + " and " + str(latePostContrastNum)
  progressBar.value = 25
  progressBar.labelText = lbltxt
  slicer.app.processEvents()
  #print("Pre-contrast image loaded to Slicer")

  #Edit 6/9/2020: Even for GE or Siemens, use Philips method of loading images into numpy array if all DCE images are in same folder

  folder1info = all_folders_info[0]
  #Loading pre-contrast image into numpy array
  #5/20/2021: Take Philips out of this if statement because Philips can have 2 folder DCE
  if (len(dce_folders) == 1):
    m,a = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,0)
  else:
    apath = exampath + "\\" + str(dce_folders[0]) 
    m,a = read_DCE_images_to_numpy.readInputToNumpy(apath)

  print("RAS to IJK Matrix")
  print(m)


  #create node for displaying pre-contrast image
  adisp = np.transpose(a,(2,1,0)) #nii needs x,y,z to have same orientation as DICOM, but for numpy array you need to return to dimension order z,y,x
  precontrast_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",prenodestr) #add this image to dropdown node with name precontrast
  slicer.util.updateVolumeFromArray(precontrast_node, adisp)
  precontrast_node.SetRASToIJKMatrix(m)

  m1 = vtk.vtkMatrix4x4()
  precontrast_node.GetIJKToRASMatrix(m1)
  print("IJK to RAS Matrix from Precontrast Node")
  print(m1)
  

  progressBar.value = 50
  progressBar.labelText = 'Pre-contrast image loaded to Slicer'
  slicer.app.processEvents()
  print("Pre-contrast image loaded to Slicer")


  #Loading early post-contrast image into numpy array
  print("-----EARLY POST-CONTRAST IMAGE-----")
  #5/20/2021: Take Philips out of this if statement because Philips can have 2 folder DCE
  if (len(dce_folders) == 1):
    m,b = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,earlyPostContrastNum)
  else:
    #Edit 11/25/2020: Add case for 2 folder DCE exams, such as those from UKCC
    if(len(dce_folders) == 2):
      m,b = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,earlyPostContrastNum)
    else:
      m,b = read_DCE_images_to_numpy.earlyOrLateImgSelect(earlyPostContrastNum,dce_folders,exampath) 

  #Edit 7/2/2020: Loading early post-contrast image directly to Slicer because want subtraction to be a button option
  #in 2nd module
  #create node for displaying early post-contrast image
  bdisp = np.transpose(b,(2,1,0)) #nii needs x,y,z to have same orientation as DICOM, but for numpy array you need to return to dimension order z,y,x
  early_post_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",earlynodestr) #add this image to dropdown node with name precontrast
  slicer.util.updateVolumeFromArray(early_post_node, bdisp)
  early_post_node.SetRASToIJKMatrix(m)

  progressBar.value = 75
  progressBar.labelText = 'Early post-contrast image loaded to Slicer'
  slicer.app.processEvents()
  print("Early post-contrast image loaded to Slicer")


  #Loading late post-contrast image into numpy array
  print("-----LATE POST-CONTRAST IMAGE-----")
  #5/20/2021: Take Philips out of this if statement because Philips can have 2 folder DCE
  if (len(dce_folders) == 1):
    m,c = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,latePostContrastNum)
  else:
    #Edit 11/25/2020: Add case for 2 folder DCE exams, such as those from UKCC
    if(len(dce_folders) == 2):
      m,c = read_DCE_images_to_numpy.readPhilipsImageToNumpy(exampath,dce_folders,fsort,latePostContrastNum)
    else:
      m,c = read_DCE_images_to_numpy.earlyOrLateImgSelect(latePostContrastNum,dce_folders,exampath) 

  #Edit 7/2/2020: Loading early post-contrast image directly to Slicer because want subtraction to be a button option
  #in 2nd module
  #create node for displaying early post-contrast image
  cdisp = np.transpose(c,(2,1,0)) #nii needs x,y,z to have same orientation as DICOM, but for numpy array you need to return to dimension order z,y,x
  late_post_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",latenodestr) #add this image to dropdown node with name precontrast
  slicer.util.updateVolumeFromArray(late_post_node, cdisp)
  late_post_node.SetRASToIJKMatrix(m)

  progressBar.value = 99
  progressBar.labelText = 'Late post-contrast image loaded to Slicer'
  slicer.app.processEvents()
  print("Late post-contrast image loaded to Slicer")

  time.sleep(1) #pause for 1 second at progress bar 99%

  progressBar.value = 100
  progressBar.labelText = 'Processing completed'
  slicer.app.processEvents()

  return precontrast_node, early_post_node, late_post_node

  

class DCE_IDandPhaseSelect(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    #Edit 7/10/2020: Change title that appears in Slicer to reflect fact that it only loads pre, early, and late images
    self.parent.title = "Module 1: Load DCE Images"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["FTV Segmentation"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["Rohan Nadkarni (UCSF Breast Imaging Research Group)"]  # TODO: replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
"""  # TODO: update with short description of the module
    self.parent.helpText += self.getDefaultModuleDocumentationLink()  # TODO: verify that the default URL is correct or change it to the actual documentation
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""  # TODO: replace with organization, grant and thanks.

#
# DCE_IDandPhaseSelectWidget
#

class DCE_IDandPhaseSelectWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    self.parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #Edit 2/11/2021: Moving exampath selection here
    self.exampath = qt.QFileDialog.getExistingDirectory(0,("Select folder for DCE-MRI exam"))


    #If using mapped drive (B: for ispy_2019 or F: for ispy2)
##    if('B:' in self.exampath or 'F:' in self.exampath):
##      #site folder name is between the 1st and 2nd slashes
##      self.sitestr = self.exampath[(int(slashinds[0])+1):int(slashinds[1])]
##      #ISPY ID folder name is between the 2nd and 3rd slashes
##      self.idstr = self.exampath[(int(slashinds[1])+1):int(slashinds[2])]
##      self.idpath = self.exampath[0:int(slashinds[2])] #full path to ispy id folder
##      #folder with visit name in it is between the 3rd and 4th slashes
##      self.visitstr = self.exampath[(int(slashinds[2])+1):int(slashinds[3])]
##
##      #studystr is ispy_2019 for B: or ispy2 for F:
##      #6/28/2021: Changing this because David said
##      #ispy_2019 is not a separate study from ispy2,
##      #it's just the name for the additional disk.
##      if('B:' in self.exampath or 'F:' in self.exampath):
##        self.studystr = 'ispy2'

    #7/6/2021: Add code to convert mapped drive letter path
    #to full path based on Andras Lasso's answer in Slicer forums.
    if(':' in self.exampath and 'C:' not in self.exampath):
      self.exampath = str(pathlib.Path(self.exampath).resolve().as_posix())
      print(self.exampath)

    #Edit 2/17/2021: Moving code to read exam details from exampath here.
    #As part of this, need to include visitstr as an input to run function
    #find indices where slashes '/' occur in exampath
    slashinds = []
    for i in range(len(self.exampath)):
      if(self.exampath[i] == '/'):
        slashinds.append(i)


    #If using full path instead of mapped drives
##    else:
    #6/28/2021: Make this code compatible with
    #directory structures that are different from
    #our MR exam directories on \\researchfiles.radiology.ucsf.edu
    if('//researchfiles' in self.exampath and ('ispy2' in self.exampath or 'ispy_2019' in self.exampath or 'acrin_6698' in self.exampath) ):
      #study folder name is between the 4th and 5th slashes
      self.studystr = self.exampath[(int(slashinds[3])+1):int(slashinds[4])]
      #Correction: ispy_2019 disk contains exams that belong to ispy2 study
      if(self.studystr == 'ispy_2019'):
        self.studystr = 'ispy2'
      #site folder name is between the 5th and 6th slashes
      self.sitestr = self.exampath[(int(slashinds[4])+1):int(slashinds[5])]
      #ISPY ID folder name is between the 6th and 7th slashes
      self.idstr = self.exampath[(int(slashinds[5])+1):int(slashinds[6])]
      self.idpath = self.exampath[0:int(slashinds[6])] #full path to ispy id folder
      #folder with visit name in it is between the 7th and 8th slashes
      self.visitstr = self.exampath[(int(slashinds[6])+1):int(slashinds[7])]
    else:
      #6/28/2021: Read info from DICOM header instead of
      #directory for 'generic' exams with other directory structures
      folders = [directory for directory in os.listdir(self.exampath) if os.path.isdir(self.exampath + "\\" + directory)]
      dcm_folder_found = 0
      for i in range(len(folders)):
        curr_path = self.exampath + "\\" + folders[i]
        curr_files = [f for f in os.listdir(curr_path) if f.endswith('.dcm')]
        curr_FILES = [f for f in os.listdir(curr_path) if f.endswith('.DCM')]
        files_noext = [f for f in os.listdir(curr_path) if f.isdigit()]

        if(len(curr_files) > 2):
          dcm1path = curr_path + "\\" + curr_files[0]
          dcm_folder_found = 1

        if(len(curr_FILES) > 2):
          dcm1path = curr_path + "\\" + curr_FILES[0]
          dcm_folder_found = 1

        if(len(files_noext) > 2):
          dcm1path = curr_path + "\\" + files_noext[0]
          dcm_folder_found = 1

        if(dcm_folder_found == 1):
          hdr_dcm1 = pydicom.dcmread(dcm1path,stop_before_pixels = True)
          try:
            self.studystr = hdr_dcm1[0x12,0x10].value #Clinical Trial Sponsor Name
          except:
            try:
              self.studystr = hdr_dcm1[0x8,0x1030].value #Study Description
            except:
              self.studystr = 'Trial Name Unknown'

          try:
            self.sitestr = hdr_dcm1[0x12,0x31].value #Clinical Trial Site Name
          except:
            try:
              self.sitestr = hdr_dcm1[0x8,0x80].value #Institution Name
            except:
              self.sitestr = 'Site Unknown'

          try:  
            self.idstr = hdr_dcm1[0x12,0x40].value #Clinical Trial Subject ID
          except:
            self.idstr = 'ID Unknown'              
            
          try:
            self.visitstr = hdr_dcm1[0x12,0x50].value #Clinical Trial Time Point ID
            #7/4/2021: Try to use directory structure if visit number not found
            #in Clinical Trial Time Point id header field
            if('1' not in self.visitstr and '2' not in self.visitstr and '3' not in self.visitstr and '4' not in self.visitstr and '5' not in self.visitstr):
              #6/29/2021: Default to 'Visit unknown',
              #change this if visit is found in exampath
              self.visitstr = 'Visit unknown'
              if('v10' in self.exampath):
                self.visitstr = 'MR1'
              
              if('v20' in self.exampath):
                self.visitstr = 'MR2'
              
              if('v30' in self.exampath):
                self.visitstr = 'MR3'
              
              if('v40' in self.exampath):
                self.visitstr = 'MR4'

          except:
            #6/29/2021: Default to 'Visit unknown',
            #change this is visit is found in exampath
            self.visitstr = 'Visit unknown'
            if('v10' in self.exampath):
              self.visitstr = 'MR1'
              
            if('v20' in self.exampath):
              self.visitstr = 'MR2'
              
            if('v30' in self.exampath):
              self.visitstr = 'MR3'
              
            if('v40' in self.exampath):
              self.visitstr = 'MR4'
              
                
    #7/26/2021: If this step fails, DICOMs in exam directory are compressed.    
    try:
      print(self.exampath)
      print(self.studystr)
      print(self.sitestr)
      print(self.idstr)
      print(self.visitstr)
    except:
      slicer.util.confirmOkCancelDisplay("Error. Please decompress all files in exam directory, then try running module again.","Compressed DICOMs Error")


    #7/16/2021: Force side-by-side axial-sagittal layout immediately after
    #processing user's MR study selection.
    slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)


    #5/18/2021: Try to retrieve study, site, visit, etc from dicom header
    #instead of exampath, because users outside UCSF will not have the
    #same directory structure
    #5/19/2021: Comment out this part when running gzipped exam
##    print("\nExam info from scan 1 DICOM header")
##    scans = [f for f in os.listdir(self.exampath) if f.isnumeric()]
##    scan1path = os.path.join(self.exampath,scans[0])
##    dcms_scan1 = os.listdir(scan1path)
##    dcm1scan1path = os.path.join(scan1path,dcms_scan1[0])
##    hdr_dcm1scan1 = pydicom.dcmread(dcm1scan1path,stop_before_pixels=True)
##    #5/18/2021: Results from UCSD 12013 MR1 suggest that I may have to
##    #find out which study based on study date instead of retrieving
##    #study directly from header field
##    studyfromhdr = hdr_dcm1scan1[0x18,0x1030].value
##    print(studyfromhdr)
##    #Group 12 is clinical trial header fields,
##    #so I'm trying to use those as much as possible.
##    sitefromhdr = hdr_dcm1scan1[0x12,0x31].value
##    print(sitefromhdr)
##    idfromhdr = hdr_dcm1scan1[0x12,0x40].value
##    print(idfromhdr)
##    visitfromhdr = hdr_dcm1scan1[0x12,0x50].value
##    print(visitfromhdr)
##    datefromhdr = hdr_dcm1scan1[0x8,0x20].value
##    print(datefromhdr)
##    print("\n")

    #7/6/2021: Add boxes to allow user to set early and late
    #post-contrast times
    
    self.earlytimelbl = qt.QLabel("\nEarly Post-Contrast Time (seconds after contrast injection):")
    self.parametersFormLayout.addRow(self.earlytimelbl)
    self.earlytime = qt.QSpinBox()
    self.earlytime.setMinimum(60)
    self.earlytime.setMaximum(180)
    self.earlytime.setSingleStep(15)
    self.earlytime.setValue(150) #set default value of 150 for earlytime spin box
    #set min and max possible values of 30 and 180 (seconds) for this spin box
    self.parametersFormLayout.addRow(self.earlytime)

    self.latetimelbl = qt.QLabel("\nLate Post-Contrast Time (seconds after contrast injection):")
    self.parametersFormLayout.addRow(self.latetimelbl)
    self.latetime = qt.QSpinBox()
    self.latetime.setMinimum(210)
    self.latetime.setMaximum(510)
    self.latetime.setSingleStep(15)
    self.latetime.setValue(450) #set default value of 450 for earlytime spin box
    #set min and max possible values of 210 and 510 (seconds) for this spin box
    self.parametersFormLayout.addRow(self.latetime)


    if('ispy' in self.exampath or 'acrin' in self.exampath):
      #position of v in visit string
      vpos = self.visitstr.find('v')
      #2/17/2021: Add exam details retrieved from exampath to DCE folder ID method selection menu.
      visnum = self.visitstr[vpos:vpos+3]
      visforlbl = 'Visit unknown'
      if(visnum == 'v10'):
        visforlbl = 'MR1'
        
      if(visnum == 'v20'):
        visforlbl = 'MR2'
        
      if(visnum == 'v25'):
        visforlbl = 'MR2.5'
        
      if(visnum == 'v30'):
        visforlbl = 'MR3'
        
      if(visnum == 'v40'):
        visforlbl = 'MR4'
        
      if(visnum == 'v50'):
        visforlbl = 'MR5'

      self.examlbl = "Choose Preferred Method of DCE Series ID for " + self.sitestr[5:] + " " + self.idstr + " " + visforlbl
    else:
      visnum = self.visitstr
      visforlbl = 'Visit unknown'

      if(visnum == 'v10'):
        visforlbl = 'MR1'
        
      if(visnum == 'v20'):
        visforlbl = 'MR2'
        
      if(visnum == 'v25'):
        visforlbl = 'MR2.5'
        
      if(visnum == 'v30'):
        visforlbl = 'MR3'
        
      if(visnum == 'v40'):
        visforlbl = 'MR4'
        
      if(visnum == 'v50'):
        visforlbl = 'MR5'

      self.examlbl = "Choose Preferred Method of DCE Series ID for " + self.sitestr + " " + self.idstr + " " + visforlbl 

    #2/11/2021: Adding interface to allow user to choose between manual and automatic folder ID
    self.DCE_ID_choose = qt.QLabel()
    self.DCE_ID_choose.setText(self.examlbl)
    self.parametersFormLayout.addRow(self.DCE_ID_choose)

    #2/11/2021: Adding checkbox for user to select automatic DCE folder ID
    self.autoIdCheckBox = qt.QCheckBox("Automatic")
    self.autoIdCheckBox.setChecked(False)
    self.parametersFormLayout.addRow(self.autoIdCheckBox)


    #2/11/2021: Adding checkbox for user to select manual DCE folder ID
    self.manualIdCheckBox = qt.QCheckBox("Manual")
    self.manualIdCheckBox.setChecked(False)
    self.parametersFormLayout.addRow(self.manualIdCheckBox)

    #2/12/2021: Adding button to submit manual DCE selections,
    #but don't actually show it unless user has selected manual option.
    self.manualSubmitButton = qt.QPushButton("Done")
    self.manualSubmitButton.setStyleSheet("margin-left:50%; margin-right:50%;") #Center align the done button for the DCE folder checkboxes
    self.manualSubmitButton.toolTip = "Submit manual DCE folder selections."
    self.manualSubmitButton.enabled = True



    #
    # Apply Button
    #
##    self.applyButton = qt.QPushButton("Run")
##    self.applyButton.toolTip = "Run the algorithm."
##    self.applyButton.enabled = False
##    self.parametersFormLayout.addRow(self.applyButton)

    # connections
##    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.autoIdCheckBox.stateChanged.connect(self.onApplyButton) #2/11/2021: connect checking the auto folder ID box to onApplyButton function
    self.manualIdCheckBox.stateChanged.connect(self.manualDCESelectMenu) #2/11/2021: connect checking the manual folder ID box to function
                                                                         #for creating (or removing) manual DCE folder selection menu.
    self.manualSubmitButton.connect('clicked(bool)',self.onApplyButton)
    
    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
##    self.onSelect()


  def cleanup(self):
    pass

##  def onSelect(self):
##    self.applyButton.enabled = True

  #2/11/2021: New function to create or remove manual DCE selection menu.
  def manualDCESelectMenu(self):
    
    if(self.manualIdCheckBox.isChecked() == True):
      #2/12/2021: If doing manual DCE folder ID, need to gunzip right now so that
      #all DICOM folders show up in the checkbox list for manual DCE folder selection.

      #Edit 7/20/2020: First thing to do is check if exam's DICOMs are gzipped
      gzipped = ident_gzipped_exam.checkForGzippedDicoms(self.exampath)

      #if exam is gzipped, gunzip it and set new exampath to gunzip folder inside of exam folder
      if(gzipped == 1):
        gzip_gunzip_pyfuncs.extractGZ(self.exampath) #create folders with unzipped DICOMs
        self.exampath = self.exampath + "\\" + "gunzipped" #set new exampath to gunzip folder inside of exam folder

      
      #First, fill header info structures for all DICOM folders in exam
      self.img_folders, self.all_folders_info = Get_header_info_all_manufacturer.fillExamFolderInfoStructures(self.exampath)

      self.dce_folders_manual = [] #Initialize array that will be filled with user's folder selections from checkboxes
      self.dce_ind_manual = [] #Indices in img_folders that correspond to user selected DCE folders
    
      #Then, add QLabel for this menu
      self.listLabel = qt.QLabel()
      self.listLabel.setText("Choose DCE folders from list below, then click 'Done' button.")
      self.listLabel.setStyleSheet("margin-left:50%; margin-right:50%;") #Center align the label for the DCE folder checkboxes
      self.parametersFormLayout.addRow(self.listLabel)

      self.listCheckBox = ["Checkbox_1", "Checkbox_2", "Checkbox_3", "Checkbox_4", "Checkbox_5",
                           "Checkbox_6", "Checkbox_7", "Checkbox_8", "Checkbox_9", "Checkbox_10",
                           "Checkbox_11", "Checkbox_12", "Checkbox_13", "Checkbox_14", "Checkbox_15",
                           "Checkbox_16", "Checkbox_17", "Checkbox_18", "Checkbox_19", "Checkbox_20",
                           "Checkbox_21", "Checkbox_22", "Checkbox_23", "Checkbox_24", "Checkbox_25",
                           "Checkbox_26", "Checkbox_27", "Checkbox_28", "Checkbox_29", "Checkbox_30",
                           "Checkbox_31", "Checkbox_32", "Checkbox_33", "Checkbox_34", "Checkbox_35",
                           "Checkbox_36", "Checkbox_37", "Checkbox_38", "Checkbox_39", "Checkbox_40",
                           "Checkbox_41", "Checkbox_42", "Checkbox_43", "Checkbox_44", "Checkbox_45",
                           "Checkbox_46", "Checkbox_47", "Checkbox_48", "Checkbox_49", "Checkbox_50"]
      allfolders_str = ''
      #Use loop to populate string that contains all DCE folder names
      for iii in range(len(self.img_folders)):
        curr_finf = self.all_folders_info[iii]
        curr_str = str(self.img_folders[iii]) + ': ' + curr_finf.serdesc
        self.listCheckBox[iii] = qt.QCheckBox(curr_str)
        self.listCheckBox[iii].setStyleSheet("margin-left:50%; margin-right:50%;") #Center align these DCE folder checkboxes
        self.listCheckBox[iii].setChecked(False)
        self.parametersFormLayout.addRow(self.listCheckBox[iii])

      #Finally, add button that will be used to submit the user's final DCE folder selections
      self.parametersFormLayout.addRow(self.manualSubmitButton)

        

    #Code for removing the manual DCE folder selection menu when Manual checkbox is unselected
    if(self.manualIdCheckBox.isChecked() == False):
      self.parametersFormLayout.removeRow(self.listLabel)

      for jjj in range(len(self.img_folders)):
        self.parametersFormLayout.removeRow(self.listCheckBox[jjj])

      self.parametersFormLayout.removeRow(self.manualSubmitButton)
        
##        allfolders_str = allfolders_str + '\n' + curr_str
##      #Once done populating this string, add it to a QLabel and display it
##      self.allDCElabel = qt.QLabel()
##      self.allDCElabel.setText(allfolders_str)
##      self.parametersFormLayout.addRow(self.allDCElabel)

          

    
  def onApplyButton(self):
    logic = DCE_IDandPhaseSelectLogic()

    if(self.autoIdCheckBox.isChecked() == True):
      self.dce_folders_manual = []
      self.dce_ind_manual = []
      logic.run(self.exampath,self.dce_folders_manual,self.dce_ind_manual,self.visitstr,self.earlytime.value,self.latetime.value)
  
    if(self.manualIdCheckBox.isChecked() == True):
      for mmm in range(len(self.img_folders)):
        if(self.listCheckBox[mmm].isChecked() == True):
          self.dce_folders_manual.append(self.img_folders[mmm])
          self.dce_ind_manual.append(mmm)
      print("Manually selected DCE folders are")
      print(self.dce_folders_manual)
      logic.run(self.exampath,self.dce_folders_manual,self.dce_ind_manual,self.visitstr,self.earlytime.value,self.latetime.value)
#
# DCE_IDandPhaseSelectLogic
#

class DCE_IDandPhaseSelectLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def run(self,exampath,dce_folders_manual,dce_ind_manual,visitstr,earlyadd,lateadd):
    """
    Run the actual algorithm
    """
    
      
      
##    if not self.isValidInputOutputData(outputVolume):
##      #slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
##      return False
    start = datetime.datetime.now()
    #Format of runstarttime is like 20200101_12:30
    runstarttime = str(start.year) + str(start.month) + str(start.day) + "_" + str(start.hour) + "_" + str(start.minute)  #Processing start time, that is included in output files or folders
    #outputfolder = "subtractions_and_MIPs_" + runstarttime

    #add exampath to parameter node
    exampath_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode","path node")
    #edit 2/16/2021: If you did manual DCE folder ID, exampath might already have "gunzipped" in it,
    #but you don't want that part of the folder path saved to the parameter node.
    #Hopefully this edit will prevent Slicer Reports from being saved inside gunzipped folder for manual
    #DCE folder ID exams with gzipped dicoms.
    if('gunzipped' in exampath):
      exampath_node.SetParameter("exampath",exampath[:-9])
    else:
      exampath_node.SetParameter("exampath",exampath)
    
    #slicer.mrmlScene.AddNode(exampath_node)


    

    #savefolder = qt.QFileDialog.getExistingDirectory(0,("Choose save location for module outputs"))
    #savepath = savefolder+"\\"+outputfolder
    
    #if not os.path.exists(savepath):
      #os.mkdir(savepath)

    #logging.info('Processing started')

    #Edit 2/12/21: Only do gunzip here if using automatic DCE folder ID.
    #If doing manual folder ID, gunzipping is already done.
    if(len(dce_folders_manual) == 0):
      #Edit 7/20/2020: First thing to do is check if exam's DICOMs are gzipped
      gzipped = ident_gzipped_exam.checkForGzippedDicoms(exampath)

      #if exam is gzipped, gunzip it and set new exampath to gunzip folder inside of exam folder
      if(gzipped == 1):
        gzip_gunzip_pyfuncs.extractGZ(exampath) #create folders with unzipped DICOMs
        exampath = exampath + "\\" + "gunzipped" #set new exampath to gunzip folder inside of exam folder

          
    #find out which visit # is the one you're processing
    #visitnum = int(visitstr[-2:]) #visit number is last 2 digits of visit folder name
    #Edit 10/30/2020: To find visitnum, need to find position of 'v' in the visit folder name
    #and take the 2 characters following that. Realized that 'v' is not necessarily at the end
    #of the visit folder name

    #6/28/2021: Make this part compatible with exams that don't use
    #\\researchfiles MR exam directory structure
    if('ispy' in exampath or 'acrin' in exampath):
      vpos = visitstr.find('v')
      visitnum = int(visitstr[vpos+1:vpos+3])
    else:
      visitnum = visitstr
    
    #Edit 7/28/2020: Must have this before you process other visits.
    #This is because by default, the 2nd module reads the 1st node's image as precontrast.
    #This means that if you add nodes for your visit after you add nodes for other visits, the wrong
    #image will be used as precontrast when computing subtraction images.
    #Call function for generating pre, early, and late nodes for the visit & exam you will do FTV processing for
    precontrast_node, early_post_node, late_post_node = loadPreEarlyLate(exampath,visitnum,1,dce_folders_manual,dce_ind_manual,earlyadd,lateadd)


  #Edit 2/17/21: Commenting out section for loading other visits because JG and TB confirmed that this is not
  #needed and it would be too difficult to allow manual folder ID for other visits.

##    visitfolders = [directory for directory in os.listdir(idpath) if os.path.isdir(idpath + "\\" + directory)]
##    print("visit folders for this exam are:")
##    print(visitfolders)
##
##    
##    #If more than one visit present in this exam, ask user if they want to load images for other visits
##    #Edit 8/7/2020: Only ask to load other visits if primary exam selected for FTV processing is NOT baseline
##    #Making this change because in Aegis, you would only load additional prior visits, not additional later visits
##    if(len(visitfolders)>1 and visitnum > 10):
##      for j in range(len(visitfolders)):
##        currvisit = visitfolders[j]
##        vpos = visitstr.find('v')
##        currvisitnum = int(currvisit[vpos+1:vpos+3])
##        #Only have option to load images to node right now if this is not the visit you are computing FTV outputs for
##        #Edit 8/7/2020: Only have option to load images to node if it's from an earlier visit than the primary exam you chose.
##        #Reason for this is same as above 8/7/2020 edit.
##        if(currvisit != visitstr and currvisitnum < visitnum):
##          if((currvisitnum % 10) == 0):
##            currvisname = 'MR' + str(int(currvisitnum/10))
##          else:
##            currvisname = 'MR' + str(float(currvisitnum/10))
##
##          yesnostr = 'Do you want to load ' + currvisname + ' images for this patient?'
##
##          #Edit 8/7/2020: Have separate yes/no dialog for each prior visit. Only load visit's images to node if user says yes
##          if(slicer.util.confirmYesNoDisplay(yesnostr)):
##            visitpath = idpath + "\\" + currvisit
##            examfolders = [directory for directory in os.listdir(visitpath) if os.path.isdir(visitpath + "\\" + directory)]
##            exampath_currvisit = visitpath + "\\" + examfolders[0]
##            #Call function for generating pre, early, and late nodes for other visits for same patient
##            #Name variables according to visit here to prevent mix-up between visits in 2nd module
##            if(currvisitnum == 10):
##              print("loading visit 10")
##              precontrast_nodev10, early_post_nodev10, late_post_nodev10 = loadPreEarlyLate(exampath_currvisit,currvisitnum,0)
##              
##            if(currvisitnum == 20):
##              print("loading visit 20")
##              precontrast_nodev20, early_post_nodev20, late_post_nodev20 = loadPreEarlyLate(exampath_currvisit,currvisitnum,0)
##              
##            if(currvisitnum == 30):
##              print("loading visit 30")
##              precontrast_nodev30, early_post_nodev30, late_post_nodev30 = loadPreEarlyLate(exampath_currvisit,currvisitnum,0)
##              
##            if(currvisitnum == 40):
##              print("loading visit 40")
##              precontrast_nodev40, early_post_nodev40, late_post_nodev40 = loadPreEarlyLate(exampath_currvisit,currvisitnum,0)
##
##            if(currvisitnum == 50):
##              print("loading visit 50")
##              precontrast_nodev50, early_post_nodev50, late_post_nodev50 = loadPreEarlyLate(exampath_currvisit,currvisitnum,0)
##
##            #Edit 7/31/2020: Just changed this one to 25 because that is the only other value I've seen
##            #Earlier I had this as else, but that made every visit (other than main one) have its images loaded to Slicer twice
##            if(currvisitnum == 25):
##              print("loading visit 25")
##              precontrast_nodev25, early_post_nodev25, late_post_nodev25 = loadPreEarlyLate(exampath_currvisit,currvisitnum,0)




    #logging.info('Processing completed')
    print("All images loaded to Slicer")
    return True


class FTVtest4Test(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_FTVtest41()

  def test_FTVtest41(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import SampleData
    SampleData.downloadFromURL(
      nodeNames='FA',
      fileNames='FA.nrrd',
      uris='http://slicer.kitware.com/midas3/download?items=5767')
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = DCE_IDandPhaseSelectLogic()
    #self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
