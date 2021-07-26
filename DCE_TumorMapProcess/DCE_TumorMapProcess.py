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
#This is the 2nd of 2 modules in the extension FTV_process_complete.
#The purpose of this module is to allow the user to
#to manually define or import ROI and omits,
#create an SER Colorized lesion segmentation,
#and generate a report.

#During ROI selection, in addition to viewing the images
#created by 1st module, the user can view MIPs or subtractions
#of these images.

import os
import unittest
import logging
import vtk, qt, ctk, slicer
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util import numpy_support
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

try:
  import numpy as np
except:
  slicer.util.pip_install('numpy')
  import numpy as np

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

import time
import datetime

try:
  import nibabel as nib
except:
  slicer.util.pip_install('nibabel')
  import nibabel as nib
  
try:
  import cv2
except:
  slicer.util.pip_install('opencv-python')
  import cv2
  
#from tr import tr
  
try:
  import nrrd
except:
  slicer.util.pip_install('pynrrd')
  import nrrd
  
import sys
import pickle
import pathlib

try:
  import yattag
except:
  slicer.util.pip_install('yattag')
  import yattag

try:
  import matplotlib
except:
  slicer.util.pip_install('matplotlib')
  import matplotlib

import matplotlib.pyplot as plt

try:
  import skimage
except:
  slicer.util.pip_install('scikit-image')
  import skimage

#imports of my functions
import compute_lps_to_rcs
import Write_to_xml
import ftv_map_gen
import ftv_map_gen_p2
import ftv_plots
import create2DimgAllFunctions
import ident_gzipped_exam
import gzip_gunzip_pyfuncs


#Function that outputs the numpy image stored in a vtkMRMLScalarVolumeNode
def getNPImgFromNode(nodename):
  prenode = slicer.util.getNode(nodename)
  img = prenode.GetImageData()
  rows,cols,slices = img.GetDimensions()
  sc = img.GetPointData().GetScalars()
  npimg = vtk_to_numpy(sc) #
  npimg = npimg.reshape(slices,rows,cols) #3/24/2020: apparently this is the correct way to reshape numpy array so that it can be viewed in slicer as NIFTI
  return npimg


#Function that returns RAS coordinate center and radius give input volume and IJK coordinate center and radius
def IJKToRASFunc(roicenter,roiradius,inputVolume):

  #Cannot directly apply transformation to roicenter, must do it this way.
  roistart = np.zeros((1,4))
  roistart[0,0:3] = roicenter-roiradius
  roistart[0,3] = 1

  roiend = np.zeros((1,4))
  roiend[0,0:3] = roicenter+roiradius
  roiend[0,3] = 1

  #Now that I'm using IJKToRAS matrix to force image orientation to be
  #the way I want, I now have to apply additional transformation to
  #roicenter and roiradius (which are now in homogeneous coordinates)
  volIJKToRASMat = vtk.vtkMatrix4x4()
  inputVolume.GetIJKToRASMatrix(volIJKToRASMat)
    
  roistart_RAS = volIJKToRASMat.MultiplyPoint(roistart[0,:])
  roistart_RAS = roistart_RAS[0:3] #remove 1 at the end
  roistart_RAS = np.array(roistart_RAS)
    
  roiend_RAS = volIJKToRASMat.MultiplyPoint(roiend[0,:])
  roiend_RAS = roiend_RAS[0:3] #remove 1 at the end
  roiend_RAS = np.array(roiend_RAS)

  roicenter_RAS = (roistart_RAS+roiend_RAS)/2

  #Edit 9/30/2020: Commenting this out because it didn't help
  #Edit 9/15/2020: Subtract 1/2 voxel size for roicenter_RAS
##  roicenter_RAS[0] = roicenter_RAS[0] - 0.5*volIJKToRASMat.GetElement(0,0)
##  roicenter_RAS[1] = roicenter_RAS[1] - 0.5*volIJKToRASMat.GetElement(1,1)
##  roicenter_RAS[2] = roicenter_RAS[2] - 0.5*volIJKToRASMat.GetElement(2,2)

  roiradius_RAS = np.abs(roiend_RAS-roicenter_RAS)

  return roicenter_RAS, roiradius_RAS



#Function that returns IJK coordinate center and radius give input volume and RAS coordinate center and radius
def RASToIJKFunc(roicenter,roiradius,inputVolume):

  volIJKToRASMat = vtk.vtkMatrix4x4()
  inputVolume.GetIJKToRASMatrix(volIJKToRASMat)

  #Edit 9/30/2020: Commenting this out because it didn't help
  #Edit 9/15/2020: Subtract 1/2 voxel size for roicenter_RAS
##  roicenter[0] = roicenter[0] - 0.5*volIJKToRASMat.GetElement(0,0)
##  roicenter[1] = roicenter[1] - 0.5*volIJKToRASMat.GetElement(1,1)
##  roicenter[2] = roicenter[2] - 0.5*volIJKToRASMat.GetElement(2,2)


  #Cannot directly apply transformation to roicenter, must do it this way.
  roistart = np.zeros((1,4))
  roistart[0,0:3] = roicenter-roiradius
  roistart[0,3] = 1

  roiend = np.zeros((1,4))
  roiend[0,0:3] = roicenter+roiradius
  roiend[0,3] = 1

  #Now that I'm using IJKToRAS matrix to force image orientation to be
  #the way I want, I now have to apply additional transformation to
  #roicenter and roiradius (which are now in homogeneous coordinates)
  volRASToIJKMat = vtk.vtkMatrix4x4()
  inputVolume.GetRASToIJKMatrix(volRASToIJKMat)
    
  roistart_IJK = volRASToIJKMat.MultiplyPoint(roistart[0,:])
  roistart_IJK = roistart_IJK[0:3] #remove 1 at the end
  roistart_IJK = np.array(roistart_IJK)
    
  roiend_IJK = volRASToIJKMat.MultiplyPoint(roiend[0,:])
  roiend_IJK = roiend_IJK[0:3] #remove 1 at the end
  roiend_IJK = np.array(roiend_IJK)

  roicenter_IJK = (roistart_IJK+roiend_IJK)/2
  roiradius_IJK = np.abs(roiend_IJK-roicenter_IJK)

  return roicenter_IJK, roiradius_IJK

##def showdialog():
##   d = qt.QDialog()
##   lbl = qt.QLabel("Choose visit to load:")
##   lbl.move(10,10)
##   b1 = qt.QPushButton("ok",d)
##   b1.move(50,50)
##   d.setWindowTitle("Choose visit to load")
##   #d.setWindowModality(Qt.ApplicationModal)
##   b1.connect('clicked(bool)',closewindow(d))
##   d.exec_()
##
##def closewindow(d):
##  print("closing window")
##  d.close()
  
#
# roiAndOmitsTest
#



class DCE_TumorMapProcess(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Module 2: Compute FTV in ROI"  # TODO: make this more human readable by adding spaces
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
# roiAndOmitsTestWidget
#


class DCE_TumorMapProcessWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """



  
  def setup(self):

##    showdialog()
    
    ScriptedLoadableModuleWidget.setup(self)





##    #Edit 8/13/2020: save initial axial position to variable
##    layoutManager = slicer.app.layoutManager()
##    red = layoutManager.sliceWidget('Red')
##    redLogic = red.sliceLogic()
##    # Save default slice offset (center axial slice of image) position to variable
##    self.axpos = redLogic.GetSliceOffset()
##
##

    #Add Axes Orientation Marker to the red and yellow slices automatically
    redNode = slicer.util.getNode('vtkMRMLSliceNodeRed')
    redNode.SetOrientationMarkerType(redNode.OrientationMarkerTypeAxes)
    
    yellowNode = slicer.util.getNode('vtkMRMLSliceNodeYellow')
    yellowNode.SetOrientationMarkerType(yellowNode.OrientationMarkerTypeAxes)

    
    

    #According to Andras, you need to disable data probe like this to prevent corner annotation from being erased
    # Disable slice annotations immediately
    slicer.modules.DataProbeInstance.infoWidget.sliceAnnotations.sliceViewAnnotationsEnabled=False
    slicer.modules.DataProbeInstance.infoWidget.sliceAnnotations.updateSliceViewFromGUI()
    # Disable slice annotations persistently (after Slicer restarts)
    settings = qt.QSettings()
    settings.setValue('DataProbe/sliceViewAnnotations.enabled', 0)


    #Edit 7/28/2020: Try reading exampath from parameter node created in 1st module
    #instead of asking user to select folder again
    #1. Retrieve node by name
    exampath_node = slicer.util.getNode("path node")
    #2. Read exampath from current node in pathSelector
    self.exampath = exampath_node.GetParameter("exampath")


    #Edit 7/29/2020: Adding exampath parsing code from createSubtractionsAndMIPs here


    #If using mapped drive (B: for ispy_2019 or F: for ispy2)
##    if('B:' in self.exampath or 'F:' in self.exampath):
##      #site folder name is between the 1st and 2nd slashes
##      self.sitestr = self.exampath[(int(slashinds[0])+1):int(slashinds[1])]
##      #ISPY ID folder name is between the 2nd and 3rd slashes
##      self.idstr = self.exampath[(int(slashinds[1])+1):int(slashinds[2])]
##      self.idpath = self.exampath[0:int(slashinds[2])] #full path to ispy id folder
##      #folder with visit name in it is between the 3rd and 4th slashes
##      visitstr = self.exampath[(int(slashinds[2])+1):int(slashinds[3])]
##
##      #studystr is ispy_2019 for B: or ispy2 for F:
##      if('B:' in self.exampath):
##        studystr = 'ispy_2019'
##      if('F:' in self.exampath):
##        studystr = 'ispy2'

    #7/6/2021: Add code to convert mapped drive letter path
    #to full path based on Andras Lasso's answer in Slicer forums.
    if(':' in self.exampath and 'C:' not in self.exampath):
      self.exampath = str(pathlib.Path(self.exampath).resolve().as_posix())
      print(self.exampath)

    #Edit 7/20/2020: First thing to do is check if exam's original DICOMs are gzipped
    self.gzipped = ident_gzipped_exam.checkForGzippedDicoms(self.exampath)
    
    #If original DICOMs are gzipped, use gunzipped folder inside exam folder
    if(self.gzipped == 1 and 'gunzipped' not in self.exampath):
      self.imagespath = self.exampath + "\\" + "gunzipped"
      #Edit 6/11/2020: Call function that does DCE folder and early/late timing identification given exampath
    #Otherwise, use exam folder
    else:
      self.imagespath = self.exampath
      #Edit 6/11/2020: Call function that does DCE folder and early/late timing identification given exampath

    #Edit 12/7/2020: Load module 1's Exam_ident_and_timing outputs from pickle workspace
    #7/9/2021: Make this path part independent of computer & directory
    #code is stored in.
    extension_path = os.path.join( os.path.dirname( __file__ ), '..' )
    wkspc_savepath = os.path.join(extension_path,'current_exam_workspace.pickle')
    with open(wkspc_savepath,'rb') as f:
      self.tempres, self.all_folders_info, self.dce_folders, self.dce_ind, self.fsort, self.studydate, self.nslice, self.earlyPostContrastNum, self.latePostContrastNum, self.earlydiffmm, self.earlydiffss, self.latediffmm, self.latediffss = pickle.load(f)

    print("dce folders")
    print(self.dce_folders)

    pre_path = os.path.join(self.imagespath,str(self.dce_folders[0]))
    pre_files = os.listdir(pre_path)
    pre_files = sorted(pre_files)
    dcm1path = os.path.join(pre_path,pre_files[0])
    hdr_dcm1 = pydicom.dcmread(dcm1path,stop_before_pixels = True)

    #7/23/2021: Need this because image reshaping different
    #when axial slices are not square (nrow != ncol).
    self.nrow_hdr = int(hdr_dcm1[0x28,0x10].value)
    self.ncol_hdr = int(hdr_dcm1[0x28,0x11].value)
    dimstr = "axial slices are " + str(self.nrow_hdr) + " x " + str(self.ncol_hdr)
    print(dimstr)



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

                            
    #Do this so you can have MR2.5 but MR1 will not be written as MR1.0
    #Edit 10/30/2020: Use new method of getting visitnum
    if('ispy' in self.exampath or 'acrin' in self.exampath):
      vpos = self.visitstr.find('v')
      visitnum = int(self.visitstr[vpos+1:vpos+3])
      #date is either the 6 digit number just before the '_' or just after the '_'
      try:
        self.datestr = self.visitstr.split('_')[0]
      except:
        self.datestr = self.visitstr.split('_')[1]

        
      if(visitnum%10 == 0):
        mrnum = int(visitnum/10)
      else:
        mrnum = float(visitnum/10)
        
      self.nodevisstr = 'MR'+str(mrnum)+' ' #v10 = MR1, v20 = MR2, etc
    else:
      self.nodevisstr = self.visitstr + ' '
      self.datestr = hdr_dcm1[0x8,0x20].value #Study Date


    #Use info from exampath to add a label at the top of the module that
    #tells which exam is being processed
    #Create the label
    #use ttlstr and ttlstr_orig because don't want multiple tumor volume values printed to red & yellow slices
    #Edit 10/30/2020: Add exam date to this string
    self.ttlstr_orig = "Study: " + self.studystr + "\nSite: " + self.sitestr + "\nPatient ID: " + self.idstr + "\nVisit: " + self.nodevisstr + "\nExam Date: " + self.datestr
    self.ttlstr = self.ttlstr_orig
    
    #Edit 10/30/2020: Easier to keep changing tumor volume printed if it's printed to the red and yellow slices.
    #Add it to current exam details after it's been computed.
    
    #Add this label to a QLabel that will appear at the top of the widget
    #Edit 10/30/2020: Exam details only written on red & yellow slice views. Tumor volume written at top of widget.
    #self.tumvolstr = "Tumor Volume Not Yet Calculated"
    #self.heading = qt.QLabel(self.tumvolstr)


    #Edit 8/13/2020: Function to edit the slice value printed to the upper right corner of image
    #Edit 9/4/2020: Print image name to bottom right of red and yellow slices too
    def updateSlicePrint(unused1=None,unused2=None): #trying 'unused' inputs to prevent error
      if(self.axchkbox.isChecked() == True):
        redprint = "Axial MIP"
        yellowprint = ""
      else:
        if(self.sagchkbox.isChecked() == True):
          redprint = ""
          yellowprint = "Sagittal MIP"
        else:
          layoutManager = slicer.app.layoutManager()

          # Save axial slice offset position to variable
          red = layoutManager.sliceWidget('Red')
          redLogic = red.sliceLogic()
          axpos = redLogic.GetSliceOffset() #update this to store new axial slice position

          #Save sagittal slice offset position to variable
          yellow = layoutManager.sliceWidget('Yellow')
          yellowLogic = yellow.sliceLogic()
          sagpos = yellowLogic.GetSliceOffset() #update this to store new sagittal slice position
      
          #Convert new slice position(s) from RAS coordinates to IJK coordinates
          pos_RAS = np.zeros((1,4)) #slice position as a homogenous RAS matrix
          pos_RAS[0,0] = sagpos #store yellow slider position
          pos_RAS[0,2] = axpos #store red slider position
          pos_RAS[0,3] = 1
          inputVolume = self.inputSelector.currentNode()
          img = inputVolume.GetImageData()
          rows,cols,slices = img.GetDimensions() #need to use cols & slices in the printed statements
          volRASToIJKMat = vtk.vtkMatrix4x4()
          inputVolume.GetRASToIJKMatrix(volRASToIJKMat)
          pos_IJK = volRASToIJKMat.MultiplyPoint(pos_RAS[0,:])

          #7/4/2021: Add case for non-bilateral images
          #based on study description

          #Update the red slice text with new axial position
          redprint = "Slice " + str(1 + round(pos_IJK[2])) + " out of " + str(slices)

          #Update the yellow slice text with new axial position
          yellowprint = "Slice " + str(1 + round(pos_IJK[0])) + " out of " + str(cols)


          #Edit 9/4/2020: Add name of current image to bottom right of red and yellow slices
          imgname = self.inputSelector.currentNode().GetName()
           
          #Update red and yellow slices with new text
          viewax = slicer.app.layoutManager().sliceWidget('Red').sliceView()
          viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,redprint)
          viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,imgname)
          viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.LowerLeft,self.ttlstr) #Edit 10/30/2020: add exam details here
          viewax.forceRender()
              
          viewsag = slicer.app.layoutManager().sliceWidget('Yellow').sliceView()
          viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,yellowprint)
          viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,imgname)
          viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.LowerLeft,self.ttlstr) #Edit 10/30/2020: add exam details here
          viewsag.forceRender()


    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)


    #Edit 8/13/2020
    #From: https://www.slicer.org/wiki/Documentation/Nightly/ScriptRepository#Get_DataProbe_text
    #You can get the mouse location in pixel coordinates along with the pixel value at the mouse
    #by hitting the '.' (period) key in a slice view after pasting in the following code. 
##    def printDataProbe():
##      infoWidget = slicer.modules.DataProbeInstance.infoWidget
##      for layer in ('B', 'F', 'L'):
##        print(infoWidget.layerNames[layer].text, infoWidget.layerIJKs[layer].text, infoWidget.layerValues[layer].text)
##
##    s = qt.QShortcut(qt.QKeySequence('.')), mainWindow())
##    s.connect('activated()', printDataProbe)    



    #Edit 5/18/2020: Make user choose exam path as soon as they load the module
    #This avoids repeated exam path selection
    #self.exampath = qt.QFileDialog.getExistingDirectory(0,("Select folder for DCE-MRI exam"))
    #alternate method of setting exampath below doesn't work--at least not yet.
    #self.exampath = exampath #Edit 7/24/2020: Directly using exampath selected in 1st module

    #parametersFormLayout.addRow(self.heading)


    #Edit 11/25/2020: Use parameter nodes to retrieve outputs of exam_ident_and_timing from 1st module
    
    #1. Retrieve node by name
##    tempres_node = slicer.util.getNode("tempres node")
##    #2. Read parameter from current node
##    tempres = tempres_node.GetParameter("tempres")
##
##    #1. Retrieve node by name
##    allfolders_node = slicer.util.getNode("all folders node")
##    #2. Read parameter from current node
##    all_folders_info = allfolders_node.GetParameter("all_folders_info")
##
##    #1. Retrieve node by name
##    dcefolders_node = slicer.util.getNode("dce folders node")
##    #2. Read parameter from current node
##    dce_folders = dcefolders_node.GetParameter("dce_folders")
##
##    #1. Retrieve node by name
##    dceind_node = slicer.util.getNode("dce indices node")
##    #2. Read parameter from current node
##    dce_ind = dceind_node.GetParameter("dce_ind")
##
##    #1. Retrieve node by name
##    fsort_node = slicer.util.getNode("fsort node")
##    #2. Read parameter from current node
##    fsort = fsort_node.GetParameter("fsort")
##
##    #1. Retrieve node by name
##    studydate_node = slicer.util.getNode("studydate node")
##    #2. Read parameter from current node
##    studydate = studydate_node.GetParameter("studydate")
##
##    #1. Retrieve node by name
##    nslice_node = slicer.util.getNode("nslice node")
##    #2. Read parameter from current node
##    nslice = nslice_node.GetParameter("nslice")
##
##    #1. Retrieve node by name
##    earlyPostContrastNum_node = slicer.util.getNode("earlyPostContrastNum node")
##    #2. Read parameter from current node
##    earlyPostContrastNum = earlyPostContrastNum_node.GetParameter("earlyPostContrastNum")
##
##    #1. Retrieve node by name
##    earlydiffmm_node = slicer.util.getNode("earlydiffmm node")
##    #2. Read parameter from current node
##    earlydiffmm = earlydiffmm_node.GetParameter("earlydiffmm")
##
##    #1. Retrieve node by name
##    earlydiffss_node = slicer.util.getNode("earlydiffss node")
##    #2. Read parameter from current node
##    earlydiffss = earlydiffss_node.GetParameter("earlydiffss")
##
##    #1. Retrieve node by name
##    latePostContrastNum_node = slicer.util.getNode("latePostContrastNum node")
##    #2. Read parameter from current node
##    latePostContrastNum = latePostContrastNum_node.GetParameter("latePostContrastNum")
##
##    #1. Retrieve node by name
##    latediffmm_node = slicer.util.getNode("latediffmm node")
##    #2. Read parameter from current node
##    latediffmm = latediffmm_node.GetParameter("latediffmm")
##
##    #1. Retrieve node by name
##    latediffss_node = slicer.util.getNode("latediffss node")
##    #2. Read parameter from current node
##    latediffss = latediffss_node.GetParameter("latediffss")

    #End edit 11/25/2020



    #7/6/2021: Add exception for DCE series with number and letters in name,
    #like Duke TCIA.
    try:
      pre_folder = int(self.dce_folders[0])
    except:
      pre_folder = self.dce_folders[0]
      
    #call function that returns affine LPS to RCS matrix and its inverse by reading 1 DICOM slice's header in pre-contrast (or entire DCE) folder
    #Edit 11/30/2020: If you have pre and multi-volume post (like UKCC), use the post folder for this function
    if(len(self.dce_folders) == 2):
      post_folder = int(self.dce_folders[1])
      self.aff_mat,self.aff_inv_mat = compute_lps_to_rcs.computeAffineAndAffineInverse(self.imagespath,post_folder,self.nslice,self.fsort)      
    else:
      self.aff_mat,self.aff_inv_mat = compute_lps_to_rcs.computeAffineAndAffineInverse(self.imagespath,pre_folder,self.nslice,self.fsort)


    print("My affine matrix")
    print(self.aff_mat)

    f1info = self.all_folders_info[0]
    self.manufacturer = f1info.manufacturer    
  

    #By default, use side by side view that only shows axial and sagittal views
    slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)

    #save datetime of run to class object for file save name purposes
    #define xml file savename with datetime in it
    start = datetime.datetime.now()
    #Format of runstarttime is like 20200101_12:30
    self.runstarttime = str(start.year) + str(start.month) + str(start.day) + "_" + str(start.hour) + "_" + str(start.minute)  #Processing start time, that is included in output files or folders

    

    #Edit 9/11/2020: Test out adding a slider (which does nothing for now)
##    self.slider1 = ctk.ctkSliderWidget()
##    self.slider1.minimum = 0
##    self.slider1.maximum = 100
##    self.slider1.value = 50
##    parametersFormLayout.addRow("Slider1: ", self.slider1)




        #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Select image to view." )
    parametersFormLayout.addRow("Displayed Image: ", self.inputSelector)






    #Create ROI and add it to scene
    #Edit 6/8/2020: Make ROI always at center of image, regardless of image size
    inputVolume = self.inputSelector.currentNode()
    img = inputVolume.GetImageData()
    rows,cols,slices = img.GetDimensions()
    roicenter = np.array([int(rows/2),int(cols/2),int(slices/2)])
    roiradius = np.array([100,100,40])

    inputVolume.GetDisplayNode().SetAutoWindowLevel(False) #Do this to prevent auto window leveling
    self.windowmin = inputVolume.GetDisplayNode().GetWindowLevelMin()
    self.windowmax = inputVolume.GetDisplayNode().GetWindowLevelMax()


    #Read current node name
    self.currentnodename = self.inputSelector.currentNode().GetName()

    #Edit 7/2/2020: automatically save numpy array for pre-contrast image
    #Edit 7/29/2020: save pre-contrast image from every visit that is loaded into a numpy array
    try:
      self.a1 = getNPImgFromNode("MR1 pre-contrast")
    except:
      print("MR1 doesn't exist or is not loaded to Slicer")
      
    try:
      self.a2 = getNPImgFromNode("MR2 pre-contrast")
    except:
      print("MR2 doesn't exist or is not loaded to Slicer")

    try:
      self.a25 = getNPImgFromNode("MR2.5 pre-contrast")
    except:
      print("MR2.5 doesn't exist or is not loaded to Slicer")
    
    try:
      self.a3 = getNPImgFromNode("MR3 pre-contrast")
    except:
      print("MR3 doesn't exist or is not loaded to Slicer")
      
    try:
      self.a4 = getNPImgFromNode("MR4 pre-contrast")
    except:
      print("MR4 doesn't exist or is not loaded to Slicer")

    try:
      self.a5 = getNPImgFromNode("MR5 pre-contrast")
    except:
      print("MR5 doesn't exist or is not loaded to Slicer")

    #7/3/2021: Add another case for when visit is not found.
    try:
      if('MR' not in self.visitstr):
        precontraststr = self.visitstr + " pre-contrast"
        self.a = getNPImgFromNode(precontraststr)
    except:
      print("pre-contrast with unknown visit doesn't exist")


    #Edit 11/6/2020: Repeat part of updateSlicePrint code here because apparently just calling the function without connecting it to slice scrolling doesn't work
    layoutManager = slicer.app.layoutManager()
    #Save axial slice offset position to variable
    red = layoutManager.sliceWidget('Red')
    redLogic = red.sliceLogic()
    axpos = redLogic.GetSliceOffset() #update this to store new axial slice position

    #Save sagittal slice offset position to variable
    yellow = layoutManager.sliceWidget('Yellow')
    yellowLogic = yellow.sliceLogic()
    sagpos = yellowLogic.GetSliceOffset() #update this to store new sagittal slice position

    #Convert new slice position(s) from RAS coordinates to IJK coordinates
    pos_RAS = np.zeros((1,4)) #slice position as a homogenous RAS matrix
    pos_RAS[0,0] = sagpos #store yellow slider position
    pos_RAS[0,2] = axpos #store red slider position
    pos_RAS[0,3] = 1
    inputVolume = self.inputSelector.currentNode()
    img = inputVolume.GetImageData()
    rows,cols,slices = img.GetDimensions() #need to use cols & slices in the printed statements
    volRASToIJKMat = vtk.vtkMatrix4x4()
    inputVolume.GetRASToIJKMatrix(volRASToIJKMat)
    pos_IJK = volRASToIJKMat.MultiplyPoint(pos_RAS[0,:])

    #Update the red slice text with new axial position
    redprint = "Slice " + str(1 + round(pos_IJK[2])) + " out of " + str(slices)

    #Update the yellow slice text with new axial position
    yellowprint = "Slice " + str(1 + round(pos_IJK[0])) + " out of " + str(cols)

    #Edit 9/4/2020: Add name of current image to bottom right of red and yellow slices
    imgname = self.inputSelector.currentNode().GetName()
     
    #Update red and yellow slices with new text
    viewax = slicer.app.layoutManager().sliceWidget('Red').sliceView()
    viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,redprint)
    viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,imgname)
    viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.LowerLeft,self.ttlstr) #Edit 10/30/2020: add exam details here
    viewax.forceRender()
        
    viewsag = slicer.app.layoutManager().sliceWidget('Yellow').sliceView()
    viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,yellowprint)
    viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,imgname)
    viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.LowerLeft,self.ttlstr) #Edit 10/30/2020: add exam details here
    viewsag.forceRender()
    #END repeat of updateSlicePrint code



    #Call function that returns RAS coordinate center and radius give input volume and IJK coordinate center and radius
    roicenter_RAS, roiradius_RAS = IJKToRASFunc(roicenter,roiradius,inputVolume)
    
    self.roi = slicer.vtkMRMLAnnotationROINode()
    self.roi.SetXYZ(roicenter_RAS)
    self.roi.SetRadiusXYZ(roiradius_RAS)
    slicer.mrmlScene.AddNode(self.roi)

    #Edit 9/11/2020: Use sliders for window level and window width

##    self.slider1 = ctk.ctkSliderWidget()
##    self.slider1.minimum = 0
##    self.slider1.maximum = 10000
##    #self.slider1.value = 50
##    parametersFormLayout.addRow("Window Width: ", self.slider1)
##
##
##    self.slider2 = ctk.ctkSliderWidget()
##    self.slider2.minimum = 0
##    self.slider2.maximum = 8000
##    #self.slider1.value = 50
##    parametersFormLayout.addRow("Window Center: ", self.slider2)




    #Edit 8/11/2020: Try QSpinBox for setting window level
    #Changing this to reflect what "window" actually means
##    self.disprangelbl = qt.QLabel("\nWindow Width:")
##    parametersFormLayout.addRow(self.disprangelbl)
##    self.disprange = qt.QSpinBox()
##    #set min and max possible values for this spin box
##    self.disprange.setMinimum(0)
##    self.disprange.setMaximum(10000)
##    parametersFormLayout.addRow(self.disprange)

    #Changing this to reflect what "level" actually means
##    self.dispctrlbl = qt.QLabel("\nWindow Center:")
##    parametersFormLayout.addRow(self.dispctrlbl)
##    self.dispctr = qt.QSpinBox()
##    #set min and max possible values for this spin box
##    self.dispctr.setMinimum(0)
##    self.dispctr.setMaximum(10000)
##    parametersFormLayout.addRow(self.dispctr)

    self.switchimage = True #variable to prevent calling of adjust window level function while changing image
                            #displayed from dropdown menu


    #Edit 9/29/2020: Read values from display node to update sliders, but do not auto window level.
    
    #Revert back to AutoWindowLevel and set values in the spin boxes accordingly
    #inputVolume.GetDisplayNode().SetAutoWindowLevel(True) #Revert back to auto window level whenever you switch images
##    self.slider1.value = ( float( inputVolume.GetDisplayNode().GetWindow() ) ) 
##    self.slider2.value = ( float( inputVolume.GetDisplayNode().GetLevel() ) ) 
##    self.dfltrange = float(self.slider1.value)
##    self.dfltctr = float(self.slider2.value)


    #Checkbox for subtraction
    self.subtractCheckBox = qt.QCheckBox("Show Subtraction")
    self.subtractCheckBox.setChecked(False)
    parametersFormLayout.addRow(self.subtractCheckBox)


    #Checkboxes for axial and sagittal MIPs

    #Axial MIP checkbox
    self.axchkbox = qt.QCheckBox("Show Axial MIP")
    self.axchkbox.setChecked(False)
    parametersFormLayout.addRow(self.axchkbox)

    #Sagittal MIP checkbox
    self.sagchkbox = qt.QCheckBox("Show Sagittal MIPs")
    self.sagchkbox.setChecked(False)
    parametersFormLayout.addRow(self.sagchkbox)

    
    #Edit 7/24/2020: Turn lesion segmenting into a checkbox option instead of a button
    #with no option to remove SER lesion segmentation after it has been overlayed on
    #image
    self.segmentLesionCheckBox = qt.QCheckBox("Show Lesion Segmentation")
    self.segmentLesionCheckBox.setChecked(False)
    parametersFormLayout.addRow(self.segmentLesionCheckBox)


    #Try to make axial & sagittal MIP checkboxes mutually exclusive
    #by adding them to a QButtonGroup
    #-- don't do this because it allows selection of one or other,
    #not both, and not neither
##    self.chkgrp = qt.QButtonGroup()
##    self.chkgrp.addButton(self.axchkbox,1)
##    self.chkgrp.addButton(self.sagchkbox,2)

    
    #Initialize empty object elements that will be used to store
    #ROI and omit region(s) centers and radii. Max of 5 omit regions.
    self.roiradius = np.zeros((1,3))
    self.roicenter = np.zeros((1,3))

    self.omitradii = np.zeros((5,3))
    self.omitcenters = np.zeros((5,3))


    # 7/21/2021: Go to ROI Center button
    #
    self.goToROICenterButton = qt.QPushButton("Go to ROI Center")
    self.goToROICenterButton.toolTip = "Go to ROI Center"
    self.goToROICenterButton.enabled = True
    parametersFormLayout.addRow(self.goToROICenterButton)


    # Import ROI Button
    #
    self.importROIButton = qt.QPushButton("Import ROI and Omits from File")
    self.importROIButton.toolTip = "Import ROI From File"
    self.importROIButton.enabled = True
    parametersFormLayout.addRow(self.importROIButton)



    # Add Omit Region Button
    #
    self.addOmitButton = qt.QPushButton("Add Omit Region")
    self.addOmitButton.toolTip = "Add Omit Region"
    self.addOmitButton.enabled = True
    parametersFormLayout.addRow(self.addOmitButton)

    #Edit 6/30/2020: Disabling the invert ROI position button because it is no longer needed


    #Invert ROI and Omits Axial Position Button -- only use if you import ROI and omits and it doesn't appear to match with tumor position
    #This means images have z-order reversed, and you have to account for this in report
##    self.imgReversed = 0 #this is how you check if axial slices are in reverse order for report labeling purposes
##    self.invertROIposButton = qt.QPushButton("Invert Axial Position of ROI and Omits")
##    self.invertROIposButton.toolTip = "Invert Axial Position of ROI and Omits"
##    self.invertROIposButton.enabled = True
##    parametersFormLayout.addRow(self.invertROIposButton)


    #Add Omit region count
    self.omitCount = 0
    
    #5 omit regions max. These are stored in distinct attributes of widget class object
    self.omit1 = slicer.vtkMRMLAnnotationROINode()
    self.omit2 = slicer.vtkMRMLAnnotationROINode()
    self.omit3 = slicer.vtkMRMLAnnotationROINode()
    self.omit4 = slicer.vtkMRMLAnnotationROINode()
    self.omit5 = slicer.vtkMRMLAnnotationROINode()


    # Save Region to File Button
    #
    self.saveRegionButton = qt.QPushButton("Save ROI and Omits to File")
    self.saveRegionButton.toolTip = "Save ROI and Omits to File"
    self.saveRegionButton.enabled = True
    parametersFormLayout.addRow(self.saveRegionButton)




    
    #Get default values of roi radius and center, because these will be used in onApplyButton for unused omits
    print("default center and radius are:")
    self.roicenter_default = [0,0,0]
    self.omit1.GetXYZ(self.roicenter_default)
    print(self.roicenter_default)
    
    self.roiradius_default = [0,0,0]
    self.omit1.GetRadiusXYZ(self.roiradius_default)
    print(self.roiradius_default)



      #
    # Apply Button
    #

    #Edit 7/10/2020: Changing button label to reflect fact that
    #based on Aegis, viewing of precontrast ROI, PE ROI, & Tumor Mask ROI
    #is probably a feature that is not needed.
    self.applyButton = qt.QPushButton("Generate Report")
    self.applyButton.toolTip = "Generate Report."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)
    


    #Edit 10/5/2020: Add spin boxes (and corresponding labels) for adjusting BKG thresh % and minconnpix and PE thresholds
    self.bkgthreshlbl = qt.QLabel("\nPercentage Background Threshold:")
    parametersFormLayout.addRow(self.bkgthreshlbl)
    self.bkgthresh = qt.QSpinBox()
    #set min and max possible values of 10 and 100 (percent) for this spin box
    self.bkgthresh.setMinimum(10)
    self.bkgthresh.setMaximum(100)
    self.bkgthresh.setSingleStep(10)
    self.bkgthresh.setValue(60) #set default value of 60 for PE threshold spin box
    parametersFormLayout.addRow(self.bkgthresh)

    self.minconnthreshlbl = qt.QLabel("\nMinimum Neighbor Count:")
    parametersFormLayout.addRow(self.minconnthreshlbl)
    self.minconnthresh = qt.QSpinBox()
    #set min and max possible values of 1 and 10 (percent) for this spin box
    self.minconnthresh.setMinimum(1)
    self.minconnthresh.setMaximum(10)
    self.minconnthresh.setValue(4) #set default value of 4 (updated 4/26/21) for min connected neighbors spin box
    parametersFormLayout.addRow(self.minconnthresh)

    self.pethreshlbl = qt.QLabel("\nPeak Enhancement Threshold:")
    parametersFormLayout.addRow(self.pethreshlbl)
    self.pethresh = qt.QSpinBox()
    #set min and max possible values of 10 and 150 for PE threshold spin box
    self.pethresh.setMinimum(10)
    self.pethresh.setMaximum(150)
    self.pethresh.setSingleStep(10)
    self.pethresh.setValue(70) #set default value of 70 for PE threshold spin box
    parametersFormLayout.addRow(self.pethresh)

    #7/8/2021: Make this part work no matter which directory & computer
    #the modules are stored in.
    #The full path to DCE_TumorMapProcess.py is
    #automatically stored in __file__ so I'm using that.
    pathtoscript,scriptname = os.path.split(os.path.realpath(__file__))
    colorbar_path = os.path.join(pathtoscript,'ser_colorbar.png')
    
    #SER Color Range image label -- this method works!
    #5/26/2020: For now, don't worry about resizing the SER Colormap
    pixmap = qt.QPixmap(colorbar_path)
    #pixmap.scaledToHeight(50) --this doesn't actually do the resize either
    #smaller_pixmap = pixmap.scaled(32, 32, qt.KeepAspectRatio, qt.FastTransformation) -- gives syntax error for KeepAspectRatio
    self.label = qt.QLabel(self)
    #self.label.resize(25,50) --didn't work either
    self.label.setPixmap(pixmap)
    #self.label.setFixedSize(150,200) --this doesn't help, just cuts off image. Need other way to resize.
    parametersFormLayout.addRow(self.label)

    #Labels for SER Color Percentage Distribution RuntimeError                         
    self.TumorPctDistLbl = qt.QLabel()
    self.TumorPctDistLbl.setText("SER Color Distribution in Lesion")
    parametersFormLayout.addRow(self.TumorPctDistLbl)

    self.YellowPctLbl = qt.QLabel()
    self.YellowPctLbl.setText( "Yellow: " )
    self.YellowValue = qt.QLabel()
    parametersFormLayout.addRow(self.YellowPctLbl,self.YellowValue)

    self.RedPctLbl = qt.QLabel()
    self.RedPctLbl.setText( "Red: " )
    self.RedValue = qt.QLabel()
    parametersFormLayout.addRow(self.RedPctLbl,self.RedValue)

    self.GreenPctLbl = qt.QLabel()
    self.GreenPctLbl.setText( "Green: " )
    self.GreenValue = qt.QLabel()
    parametersFormLayout.addRow(self.GreenPctLbl,self.GreenValue)

    self.PurplePctLbl = qt.QLabel()
    self.PurplePctLbl.setText( "Purple: " )
    self.PurpleValue = qt.QLabel()
    parametersFormLayout.addRow(self.PurplePctLbl,self.PurpleValue)

    self.BluePctLbl = qt.QLabel()
    self.BluePctLbl.setText( "Blue: " )
    self.BlueValue = qt.QLabel()
    parametersFormLayout.addRow(self.BluePctLbl,self.BlueValue)

    self.TotVolLbl = qt.QLabel()
    self.TotVolLbl.setText("Total: ")
    self.TotVolValue = qt.QLabel()
    parametersFormLayout.addRow(self.TotVolLbl,self.TotVolValue)


    # connections
    #Edit 9/29/2020: try adding code to update window and level
    if(self.windowmin != inputVolume.GetDisplayNode().GetWindowLevelMin() or self.windowmax != inputVolume.GetDisplayNode().GetWindowLevelMax()):
      self.updateWindowSettingVariables
    
    #Edit: Cannot display axial & sagittal MIPs simultaneously because setting SlabNumberOfSlices for axial MIP
    #messes up sagittal MIP, and vice versa
    self.subtractCheckBox.stateChanged.connect(self.showSubtractionFromNode)
    self.axchkbox.stateChanged.connect(self.showAxialMIPFromNode)
    self.sagchkbox.stateChanged.connect(self.showSagittalMIPFromNode)
    self.segmentLesionCheckBox.stateChanged.connect(self.segmentLesion) #Edit 7/24/2020: This is now a checkbox option
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.addOmitButton.connect('clicked(bool)',self.addOmitRegion)
    self.importROIButton.connect('clicked(bool)',self.importROIFromFile)
    self.saveRegionButton.connect('clicked(bool)',self.saveRegionToXML)
    self.goToROICenterButton.connect('clicked(bool)',self.goToROICenter)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
##    self.slider1.connect('valueChanged(double)', self.adjustWindowLevel)
##    self.slider2.connect('valueChanged(double)', self.adjustWindowLevel)

##    self.dispctr.valueChanged.connect(self.adjustWindowLevel)
##    self.disprange.valueChanged.connect(self.adjustWindowLevel)

    
    #Call function to update text giving axial slice position whenever red slice position is adjusted
    slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed').AddObserver(vtk.vtkCommand.ModifiedEvent,updateSlicePrint)
    #Call function to update text giving sagittal slice position whenever yellow slice position is adjusted
    slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeYellow').AddObserver(vtk.vtkCommand.ModifiedEvent,updateSlicePrint)

      
##    self.invertROIposButton.connect('clicked(bool)',self.invertROIpos)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()


    

  def cleanup(self):
    pass


  def updateWindowSettingVariables(self):
    print("updating variables to match display setttings")
    self.windowmin = inputVolume.GetDisplayNode().GetWindowLevelMin()
    self.windowmax = inputVolume.GetDisplayNode().GetWindowLevelMax()
    

  #Edit 8/11/2020: Function for adjusting window level according to user settings from spin boxes
##  def adjustWindowLevel(self):
##    #only adjust window level if changing of window or level is not during dropdown menu image switch 
##    if(self.switchimage == False):
##      inputVolume = self.inputSelector.currentNode()
##      window = float(self.slider1.value)
##      level = float(self.slider2.value)
##      inputVolume.GetDisplayNode().SetAutoWindowLevel(False) #Disable Auto window level
##      inputVolume.GetDisplayNode().SetWindowLevel(window,level)


  def makeSubImg(self):
    #read in numpy array from node and compute subtraction image
    inputVolume = self.inputSelector.currentNode()
    img = inputVolume.GetImageData()
    rows,cols,slices = img.GetDimensions()
    sc = img.GetPointData().GetScalars()
    postimg = vtk_to_numpy(sc)
    if(self.nrow_hdr == self.ncol_hdr):
      postimg = postimg.reshape(slices,rows,cols)
    else:
      postimg = postimg.reshape(slices,cols,rows)

    #1st set of if statements: Choose the correct pre-contrast image
    #to subtract based on which visit's image you are viewing
    if('MR1' in self.currentnodename):
      a = self.a1
    if('MR2 ' in self.currentnodename): #Edit 8/14/2020: add space to end of MR2 so that MR2.5 doesn't read this if and get an error
      a = self.a2
    if('MR2.5' in self.currentnodename):
      a = self.a25
    if('MR3' in self.currentnodename):
      a = self.a3
    if('MR4' in self.currentnodename):
      a = self.a4
    if('MR5' in self.currentnodename):
      a = self.a5
    #7/3/2021: Add another case for when visit is unknown
    if('MR' not in self.currentnodename):
      a = self.a

    #2nd set of if statements: Compute subtraction image,
    #subtraction axial MIP, or subtraction sagittal MIP
    if(self.axchkbox.isChecked() == True):
        
      #if input is axial MIP, compute pre-contrast axial MIP, then use it to make a subtraction axial MIP

      if(self.nrow_hdr != self.ncol_hdr):
        a = a.reshape(slices,cols,rows)

      #1. compute pre-contrast axial MIP
      a = np.transpose(a,(2,1,0)) #change from z,y,x to x,y,z
      pre_ax_mip_1slc = create2DimgAllFunctions.makeMIP(a)
      pre_ax_mip = np.zeros((a.shape))
        
      #project across all slices
      for zval in range(a.shape[2]):
        pre_ax_mip[:,:,zval] = pre_ax_mip_1slc

      pre_ax_mip = np.transpose(pre_ax_mip,(2,1,0)) #change from x,y,z back to z,y,x

      #2. compute subtraction axial MIP
      subimg = postimg-pre_ax_mip

      #3. name node
      subnodename = self.currentnodename[0:-9] + "subtraction " + self.currentnodename[-9:]
      #7/22/2021: See if this line helps
      self.axmipnodename = subnodename

    else:
      if(self.sagchkbox.isChecked() == True):
        #if input is sagittal MIP, compute pre-contrast sagittal MIP, then use it to make a subtraction sagittal MIP

        if(self.nrow_hdr != self.ncol_hdr):
          a = a.reshape(slices,cols,rows)

        #1. compute pre-contrast sagittal MIP
        a = np.transpose(a,(2,1,0)) #change from z,y,x to x,y,z
        mid_slc = int(a.shape[0]/2)
        atrans = np.transpose(a,(1,2,0)) #change from x,y,z to y,z,x
          
        pre_sag_mip = np.zeros((a.shape)) #MIP image has 1 MIP for left half of image and another MIP for right half of image
        pre_sag_mip_1slc_h1 = create2DimgAllFunctions.makeMIP(atrans[:,:,0:mid_slc]) #x is last index

        #project across first half of x-slices (the rest are black)
        for xval in range(0,mid_slc):
          pre_sag_mip[xval,:,:] = pre_sag_mip_1slc_h1

        pre_sag_mip_1slc_h2 = create2DimgAllFunctions.makeMIP(atrans[:,:,mid_slc:]) #x is last index
        #project across second half of x-slices
        for xval in range(mid_slc,a.shape[0]):
          pre_sag_mip[xval,:,:] = pre_sag_mip_1slc_h2

        pre_sag_mip = np.transpose(pre_sag_mip,(2,1,0)) #change from x,y,z to z,y,x

        #2. compute subtraction sagittal MIP
        subimg = postimg-pre_sag_mip

        #3. name node
        subnodename = self.currentnodename[0:-12] + "subtraction " + self.currentnodename[-12:]
        #7/22/2021: See if this line helps
        self.sagmipnodename = subnodename
        
      else:
        if(self.nrow_hdr != self.ncol_hdr):
          a = a.reshape(slices,cols,rows)
          
        #If input is not a MIP, just subtract pre-contrast original image from input to display what you want
        subimg = postimg-a
          
        #name node
        subnodename = self.currentnodename + " subtraction"
        
    #Edit 7/23/2020: Set negative values in subtraction image to 0
    subimg = np.where(subimg>0,subimg,0)
          
    return inputVolume, subimg, subnodename


  def showSubtractionFromNode(self):

    #if subtraction box is checked, the subtraction image must be created
    if(self.subtractCheckBox.isChecked() == True):
      self.newSubSelection = 0

      print("Showing post-contrast subtraction image")
      self.inputToSubtract = self.inputSelector.currentNode() #use this to return to input when you uncheck box
      self.currentnodename = self.inputSelector.currentNode().GetName() #use this to check if input is early or late, and if input is MIP or not

      #call function to choose appropriate image and node name for subtraction display
      inputVolume, subimg, subnodename = self.makeSubImg()

      #display new subtraction image
      self.subimg_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",subnodename) #add this image to dropdown node
      slicer.util.updateVolumeFromArray(self.subimg_node, subimg)
      #add transformation matrix to node
      m = vtk.vtkMatrix4x4()
      inputVolume.GetRASToIJKMatrix(m)
      self.subimg_node.SetRASToIJKMatrix(m)
      #display subtraction image node
      self.inputSelector.setCurrentNode(self.subimg_node)

      if("axial MIP" in self.currentnodename):
        #if input was axial MIP, you have to manual re-check axial MIP box because changing image unchecks it
        self.axchkbox.setChecked(True)
        
      if("sagittal MIP" in self.currentnodename):
        #if input was sagittal MIP, you have to manual re-check sagittal MIP box because changing image unchecks it
        self.sagchkbox.setChecked(True)

    #If box is unchecked, revert back to showing input image
    else:
      print("Showing original post-contrast image")
      #Switch back to input image to MIP when appropriate
      #Edit 10/1/2020: Try Csaba Pinter's method for displaying image selected from dropdown menu or checkbox
      #I tried it and it does in fact allow you to switch images while retaining the zoom and slice position chosen
      #by the user in the prior image
      sliceCompositeNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
      for sliceCompositeNode in sliceCompositeNodes:
        sliceCompositeNode.SetBackgroundVolumeID(self.inputToSubtract.GetID())

    
      slicer.mrmlScene.RemoveNode(self.subimg_node)

      self.inputSelector.setCurrentNode(self.inputToSubtract)
      
    self.startInputVolume = self.inputSelector.currentNode() #volume to use for window leveling
    self.subwindow = self.startInputVolume.GetDisplayNode().GetWindow()
    self.sublevel = 0.5*(self.startInputVolume.GetDisplayNode().GetWindowLevelMin() + self.startInputVolume.GetDisplayNode().GetWindowLevelMax())



  #Function for finding axial slice with most SER values > 1.1 and <= 3
  #view: value is 'ax' or 'sag'
  def findMaxSERSlice(self,view):
    if(view == 'ax'):
      ser_test = self.ser #if axial, use ser image in its original orientation
    if(view == 'sag'):
      ser_test = np.copy(self.ser)
      ser_test = np.transpose(ser_test,(2,1,0)) #if sagittal transpose ser image from x,y,z to z,y,x
        
    slc_maxSER = 0
    maxhighSERs = 0
    #iterate through all slices, and update each time you reach a new max of most high SER values
    for i in range(ser_test.shape[2]):
      ser_slc = ser_test[:,:,i]
      ser_lbmask = (ser_slc > 1.1)
      ser_ubmask = (ser_slc <= 3)
      ser_combmask = ser_lbmask*ser_ubmask
      highSERs = np.sum(ser_combmask)

      if(highSERs > maxhighSERs):
        maxhighSERs = highSERs
        slc_maxSER = i

    return slc_maxSER


  def setAxialSliderPositionToROICenter(self,inputVolume,img):
  #function to automatically set MIP axial slice position to center axial slice of ROI



      #Edit 7/8/2020: Add code to set automatically set MIP axial slice position to center axial slice of ROI
      #Edit 7/30/2020: if there is an SER color node, use the slice with most green, red, and yellow SER voxels

      #Edit 7/30/2020: Disabling this if else for now because it is giving me problems with sagittal slice position.
      #I'm sticking with just showing center slice of ROI for SER color overlayed on MIP

##      if(hasattr(self,'segment_node') == 1):
##        roizcenterIJK = self.findMaxSERSlice('ax') #Keep same name for slice for convenience
##      #otherwise, use slice at center of ROI
##      else:
        #read ROI and omit region info
    self.roi.GetXYZ(self.roicenter[0,:])
    
    self.roi.GetRadiusXYZ(self.roiradius[0,:])
    self.roiradius = np.absolute(self.roiradius)

    #convert ROI coordinates from RAS to IJK
    self.roicenter[0,:], self.roiradius[0,:] = RASToIJKFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)
    roizcenterIJK = self.roicenter[0,2]
    roizcenterIJK = int(roizcenterIJK) #center z voxel value of ROI
    print("roi z center IJK")
    print(roizcenterIJK)

    layoutManager = slicer.app.layoutManager()
    red = layoutManager.sliceWidget('Red')
    redLogic = red.sliceLogic()
    # Save default slice offset (center axial slice of image) position to variable
    imgzcenter = redLogic.GetSliceOffset()
    print("image z center RAS")
    print(imgzcenter)

    #Get spacing between adjacent slices
    volIJKToRASMat = vtk.vtkMatrix4x4()
    inputVolume.GetIJKToRASMatrix(volIJKToRASMat)
    zspacing = volIJKToRASMat.GetElement(2,2)
    print("z spacing RAS")
    print(zspacing)

    print("image z center IJK")
    print(int(img.shape[0]/2))
      
    #compute S coordinate of center axial slice of ROI using slice offset value above, spacing between adjacent slices,
    #and difference between center z slice of ROI and center z slice of image
    roizcenter = imgzcenter + zspacing*( roizcenterIJK - int(img.shape[0]/2) )
    print("ROI z center RAS")
    print(roizcenter)
    # Change slice position to center axial slice of ROI
    redLogic.SetSliceOffset(roizcenter)


  #Edit 7/8/2020: Add code to set automatically set MIP sagittal slice position to center sagittal slice of ROI
  def setSagittalSliderPositionToROICenter(self,inputVolume,img):


      #Edit 7/30/2020: if there is an SER color node, use the slice with most green, red, and yellow SER voxels

      #Edit 7/30/2020: Disabling this if else for now because it is giving me problems with sagittal slice position.
      #I'm sticking with just showing center slice of ROI for SER color overlayed on MIP
##      if(hasattr(self,'segment_node') == 1):
##        roixcenterIJK = self.findMaxSERSlice('sag') #Keep same name for slice for convenience
##      #otherwise, use slice at center of ROI
##      else:
        #read ROI info
      self.roi.GetXYZ(self.roicenter[0,:])
    
      self.roi.GetRadiusXYZ(self.roiradius[0,:])
      self.roiradius = np.absolute(self.roiradius)

      #convert ROI coordinates from RAS to IJK
      self.roicenter[0,:], self.roiradius[0,:] = RASToIJKFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)
      roixcenterIJK = self.roicenter[0,0]
      roixcenterIJK = int(roixcenterIJK) #center x voxel value of ROI
      print("ROI x center IJK")
      print(roixcenterIJK)
      
      layoutManager = slicer.app.layoutManager()
      yellow = layoutManager.sliceWidget('Yellow')
      yellowLogic = yellow.sliceLogic()
      # Save default slice offset (center sagittal slice of image) position to variable
      imgxcenter = yellowLogic.GetSliceOffset()
      print("image x center RAS")
      print(imgxcenter)

      #Get spacing between adjacent columns
      volIJKToRASMat = vtk.vtkMatrix4x4()
      inputVolume.GetIJKToRASMatrix(volIJKToRASMat)
      xspacing = volIJKToRASMat.GetElement(0,0)
      print("xspacing RAS")
      print(xspacing)

      print("image x center IJK")
      print(img.shape[2]/2)

      #compute R coordinate of center sagittal slice of ROI using slice offset value above, spacing between adjacent columns,
      #and difference between center z slice of ROI and center z slice of image
      roixcenter = imgxcenter + xspacing*( roixcenterIJK - int(img.shape[2]/2) )
      print("ROI x center RAS")
      print(roixcenter)
      # Change slice position to center sagittal slice of ROI
      yellowLogic.SetSliceOffset(roixcenter)

  #7/21/2021: New function to fix glitches in reverting back to ROI center
  def goToROICenter(self):
    #convert input volume to numpy array
    inputVolume = self.inputSelector.currentNode()
    img = inputVolume.GetImageData()
    rows,cols,slices = img.GetDimensions()
    #read inputVolume to numpy array
    sc = img.GetPointData().GetScalars()
    img_np = vtk_to_numpy(sc)
    #I tried other dimension orders but they gave "Scrambled" image
    img_np = img_np.reshape(slices,rows,cols) #3/24/2020: apparently this is the correct way to reshape numpy array so that it can be viewed in slicer as NIFTI
    slicer.util.resetSliceViews() #want to always reset slider positions to center of image after going from MIP to non-MIP

    #Edit 7/22/2021: Removed the reset slice views from the setAx and setSag functions
    #because this negates the result of the function you call first. Pasting the reset
    #here instead.
    
    #Return to side by side view that shows axial and sagittal views
    slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)

    slicer.util.resetSliceViews() #This reset before setting slice position based on ROI may ensure that MIP fills as much of screen as possible      
    self.setSagittalSliderPositionToROICenter(inputVolume,img_np) #set sagittal position to center of ROI
    self.setAxialSliderPositionToROICenter(inputVolume,img_np) #set axial position to center of ROI


  #Function for creating and showing axial MIP
  #by making numpy array and adding it to vtkMRMLScalarVolumeNode
  def showAxialMIPFromNode(self):
    #made ax_mip_node part of self object so that it isn't removed from memory every time function is finished running
    #that way, it can be called again to remove node.

    #If axial MIP box checked, create and display axial MIP
    if(self.axchkbox.isChecked() == True):
      self.currentnodename = self.inputSelector.currentNode().GetName()
      
      print("checking axial MIP box")
      self.newMIPSelection = 0

      #First, if sagittal MIP box is checked, uncheck it
      if (self.sagchkbox.isChecked() == True):
        self.sagchkbox.setChecked(False)
        #self.axchkbox.setChecked(True) #Need to do this again to override effect in onSelect, which sets both to False
      
      #if exists, remove previous copies of axial MIP node
      try:
        slicer.mrmlScene.RemoveNode(self.ax_mip_node)
        print("deleting old axial MIP")
      except:
        print("no old axial MIP to be removed")

      print("Creating and showing new axial MIP")
    
      #convert input volume to numpy array
      inputVolume = self.inputSelector.currentNode()
      self.inputToMIP = inputVolume #need to save inputVolume ID like this so it can be used next time you uncheck box
      img = inputVolume.GetImageData()
      rows,cols,slices = img.GetDimensions()
      #read inputVolume to numpy array
      sc = img.GetPointData().GetScalars()
      img_np = vtk_to_numpy(sc)
      
      #I tried other dimension orders but they gave "Scrambled" image
      if(self.nrow_hdr == self.ncol_hdr):
        img_np = img_np.reshape(slices,rows,cols) #3/24/2020: apparently this is the correct way to reshape numpy array so that it can be viewed in slicer as NIFTI
      else:
        img_np = img_np.reshape(slices,cols,rows)
        
      img_np = img_np.transpose(2,1,0) #transpose to cols,rows,slices to give same orientation as input. Can't do this earlier because image gets "scrambled"
      img_np = img_np.astype('float64')

      #axial MIP
      ax_mip_1slc = create2DimgAllFunctions.makeMIP(img_np)
      ax_mip = np.zeros((img_np.shape))

      #project across all slices
      for zval in range(img_np.shape[2]):
        ax_mip[:,:,zval] = ax_mip_1slc
        
      #Edit 6/29/2020: Make node name depend on input
      self.axmipnodename = self.currentnodename + " axial MIP"

      #create node for displaying  axial MIP        
      self.ax_mip_disp = np.transpose(ax_mip,(2,1,0)) #nii needs x,y,z to have same orientation as DICOM, but for numpy array you need to return to dimension order z,y,x
      self.ax_mip_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",self.axmipnodename) #add this image to dropdown node with name MIPaxial_earlypost
      slicer.util.updateVolumeFromArray(self.ax_mip_node, self.ax_mip_disp)
      #add transformation matrix to node
      m = vtk.vtkMatrix4x4()
      inputVolume.GetRASToIJKMatrix(m)
      self.ax_mip_node.SetRASToIJKMatrix(m)

      #make the new axial MIP the currently selected node automatically
      #self.inputSelector.setCurrentNode(self.ax_mip_node)

      #Edit 10/1/2020: Try Csaba Pinter's method for displaying image selected from dropdown menu or checkbox
      #I tried it and it does in fact allow you to switch images while retaining the zoom and slice position chosen
      #by the user in the prior image
##      sliceCompositeNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
##      for sliceCompositeNode in sliceCompositeNodes:
##        sliceCompositeNode.SetBackgroundVolumeID(self.ax_mip_node.GetID())

      #Edit 10/15/2020: Do this and let onSelect handle switching image. This makes the dropdown menu and top left corner of
      #red & yellow slices show the correct name ... Axial MIP
      self.inputSelector.setCurrentNode(self.ax_mip_node)

      #show only red slice
      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)

      #Manually set slice position to "Axial MIP" to fix glitch
      redprint = "Axial MIP"
      viewax = slicer.app.layoutManager().sliceWidget('Red').sliceView()
      viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,redprint)
      viewax.forceRender()


      #call function to set red slice position to center axial slice of ROI
      slicer.util.resetSliceViews() #This reset before setting slice position based on ROI may ensure that MIP fills as much of screen as possible      
      #7/21/2021: re-enabling this to see what happens
      self.setAxialSliderPositionToROICenter(inputVolume,self.ax_mip_disp)
      self.setSagittalSliderPositionToROICenter(inputVolume,self.ax_mip_disp)

      #Edit 7/23/2020: Only disable slice scrolling if lesion segmentation is absent
      #Edit 7/23/2020: Disabling this if statement because it causes a glitch where you are scrolling through
      #another breast image behind the MIP. This is in spite of the fact that the SER image is
      #white in non SER colored voxels

      #if(hasattr(self,'segment_node') == 0):

      #Edit 10/15/2020: Disable all references to interactorStyle.SetActionEnabled
      #because calling this seems to disable ability to use Adjust Window Level tool
      
      #Disable slice browsing with mouse and remove slider for red slice
      #Disable mouse browsing
      #interactorStyle = slicer.app.layoutManager().sliceWidget('Red').sliceView().sliceViewInteractorStyle()
      #interactorStyle.SetActionEnabled(interactorStyle.BrowseSlice, False)
      #Remove red slice slider from view
      lm = slicer.app.layoutManager()
      lm.sliceWidget('Red').sliceController().setVisible(False)


    #If axial MIP box unchecked, remove node for axial MIP, set input image as current node
    #onSelect will take care of the rest
    else:
      #Edit 6/29/2020: Setting current node to input image before removing MIP node
      #ensures that correct image is displayed after unchecking box.
      #(Before: late_post_subtraction always displayed after unchecking box)
      print("unchecking axial MIP box")
      #Switch back to input image to MIP when appropriate
      if self.newMIPSelection == 0:
        #Edit 10/1/2020: Try Csaba Pinter's method for displaying image selected from dropdown menu or checkbox
        #I tried it and it does in fact allow you to switch images while retaining the zoom and slice position chosen
        #by the user in the prior image
        sliceCompositeNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
        for sliceCompositeNode in sliceCompositeNodes:
          sliceCompositeNode.SetBackgroundVolumeID(self.inputToMIP.GetID())

        #Return to side by side view that shows axial and sagittal views
        slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)

      slicer.mrmlScene.RemoveNode(self.ax_mip_node)

      self.inputSelector.setCurrentNode(self.inputToMIP)
      self.goToROICenter

      #Enable slice browsing with mouse and display slider for red slice
      #Enable mouse browsing
      #interactorStyle = slicer.app.layoutManager().sliceWidget('Red').sliceView().sliceViewInteractorStyle()
      #interactorStyle.SetActionEnabled(interactorStyle.BrowseSlice, True)
      #Add red slice slider to view
      lm = slicer.app.layoutManager()
      lm.sliceWidget('Red').sliceController().setVisible(True)

      
      #slicer.util.resetSliceViews() #want to always reset slider positions to center of image after going from MIP to non-MIP
      #interactorStyle.SetActionEnabled(interactorStyle.BrowseSlice, True) #re-enable slice browsing
      #7/21/2021: re-enabling this to see what happens
##      self.setSagittalSliderPositionToROICenter(self.inputSelector.currentNode(),self.ax_mip_disp) #set sagittal position to center of ROI




  #Function for creating and showing axial MIP
  #by making numpy array and adding it to vtkMRMLScalarVolumeNode
  def showSagittalMIPFromNode(self):
    #made ax_mip_node part of self object so that it isn't removed from memory every time function is finished running
    #that way, it can be called again to remove node.

    

    #If sagittal MIP box checked, create and display sagittal MIP
    if(self.sagchkbox.isChecked() == True):
      self.currentnodename = self.inputSelector.currentNode().GetName()

      print("checking sagittal MIP box")
      self.newMIPSelection = 0

      #First, if axial MIP box is checked, uncheck it
      if (self.axchkbox.isChecked() == True):
        self.axchkbox.setChecked(False)
        #self.sagchkbox.setChecked(True) #Need to do this again to override effect in onSelect, which sets both to False

      
      #if exists, remove previous copies of axial MIP node
      try:
        slicer.mrmlScene.RemoveNode(self.sag_mip_node)
        print("deleting old sagittal MIP")
      except:
        print("no old sagittal MIP to be removed")

      print("Creating and showing new sagittal MIP")
    
      #convert input volume to numpy array
      inputVolume = self.inputSelector.currentNode()
      self.inputToMIP = inputVolume #need to save inputVolume ID like this so it can be used next time you uncheck box
      img = inputVolume.GetImageData()
      rows,cols,slices = img.GetDimensions()
      #read inputVolume to numpy array
      sc = img.GetPointData().GetScalars()
      img_np = vtk_to_numpy(sc)
      
      #I tried other dimension orders but they gave "Scrambled" image
      if(self.nrow_hdr == self.ncol_hdr):
        img_np = img_np.reshape(slices,rows,cols) #3/24/2020: apparently this is the correct way to reshape numpy array so that it can be viewed in slicer as NIFTI
      else:
        img_np = img_np.reshape(slices,cols,rows)

      img_np = img_np.transpose(2,1,0) #transpose to cols,rows,slices to give same orientation as input. Can't do this earlier because image gets "scrambled"
      img_np = img_np.astype('float64')
      img_np_trans = np.transpose(img_np,(1,2,0)) #y,z,x

      mid_slc = int(img_np.shape[0]/2) #center x index of original image

      sag_mip = np.zeros((img_np.shape)) #MIP image has 1 MIP for left half of image and another MIP for right half of image
      sag_mip_1slc_h1 = create2DimgAllFunctions.makeMIP(img_np_trans[:,:,0:mid_slc]) #x is last index

      #project across first half of x-slices (the rest are black)
      for xval in range(0,mid_slc):
        sag_mip[xval,:,:] = sag_mip_1slc_h1

      sag_mip_1slc_h2 = create2DimgAllFunctions.makeMIP(img_np_trans[:,:,mid_slc:]) #x is last index
      #project across second half of x-slices
      for xval in range(mid_slc,img_np.shape[0]):
        sag_mip[xval,:,:] = sag_mip_1slc_h2
        
      #Edit 6/29/2020: Make node name depend on input
      self.sagmipnodename = self.currentnodename + " sagittal MIP"

      #create node for displaying early post-contrast subtraction sagittal MIP
      self.sag_mip_disp = np.transpose(sag_mip,(2,1,0)) #nii needs x,y,z to have same orientation as DICOM, but for numpy array you need to return to dimension order z,y,x
      self.sag_mip_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",self.sagmipnodename) #add this image to dropdown node with name MIPsag_earlypost_x0tocenter
      slicer.util.updateVolumeFromArray(self.sag_mip_node, self.sag_mip_disp)
      #add transformation matrix to node
      m = vtk.vtkMatrix4x4()
      inputVolume.GetRASToIJKMatrix(m)
      self.sag_mip_node.SetRASToIJKMatrix(m)
      
      #make the new axial MIP the currently selected node automatically
      #self.inputSelector.setCurrentNode(self.sag_mip_node)


      #make the new sagittal MIP the currently selected node automatically
      #self.inputSelector.setCurrentNode(self.sag_mip_node)
      #Edit 10/1/2020: Try Csaba Pinter's method for displaying image selected from dropdown menu or checkbox
      #I tried it and it does in fact allow you to switch images while retaining the zoom and slice position chosen
      #by the user in the prior image
##      sliceCompositeNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
##      for sliceCompositeNode in sliceCompositeNodes:
##        sliceCompositeNode.SetBackgroundVolumeID(self.sag_mip_node.GetID())

      #Edit 10/15/2020: Use this method to make the new axial MIP the currently selected node automatically
      #so that dropdown menu and top left corner of red & yellow slices show the name ... sagittal MIP
      self.inputSelector.setCurrentNode(self.sag_mip_node)

      #show only yellow slice
      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView)

      #Manually set slice position to "Sagittal MIP" to fix glitch
      yellowprint = "Sagittal MIP"
      viewsag = slicer.app.layoutManager().sliceWidget('Yellow').sliceView()
      viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,yellowprint)
      viewsag.forceRender()



      #reset slice views so that you are shifting slice with center sagittal slice as starting point
      #slicer.util.resetSliceViews() #This reset before setting slice position based on ROI may ensure that MIP fills as much of screen as possible      

      slicer.util.resetSliceViews() #This reset before setting slice position based on ROI may ensure that MIP fills as much of screen as possible      
      #7/21/2021: re-enabling this to see what happens
      self.setSagittalSliderPositionToROICenter(inputVolume,self.sag_mip_disp)
      self.setAxialSliderPositionToROICenter(inputVolume,self.sag_mip_disp)
      #7/21/2021: call axial version of this too so that axial slice position after
      #unchecking sagittal MIP box is helpful for user.
##      self.setAxialSliderPositionToROICenter(inputVolume,self.ax_mip_disp)


      #Edit 7/23/2020: Only disable slice scrolling if lesion segmentation is absent
      #Edit 7/23/2020: Disabling this if statement because it causes a glitch where you are scrolling through
      #another breast image behind the MIP. This is in spite of the fact that the SER image is
      #white in non SER colored voxels
      
      #if(hasattr(self,'segment_node') == 0):
      
      #Disable slice browsing with mouse and remove slider for yellow slice
      #Disable mouse browsing
      #interactorStyle = slicer.app.layoutManager().sliceWidget('Yellow').sliceView().sliceViewInteractorStyle()
      #interactorStyle.SetActionEnabled(interactorStyle.BrowseSlice, False)
      #Remove yellow slice slider from view
      lm = slicer.app.layoutManager()
      lm.sliceWidget('Yellow').sliceController().setVisible(False)
      

    #If sagittal MIP box unchecked, remove node for sagittal MIP, set input image as current node
    #onSelect will take care of the rest
    else:
      #Edit 6/29/2020: Setting current node to input image before removing MIP node
      #ensures that correct image is displayed after unchecking box.
      #(Before: late_post_subtraction always displayed after unchecking box)
      print("unchecking sagittal MIP box")
      #Switch back to input image to MIP when appropriate
      if self.newMIPSelection == 0:
        #Edit 10/1/2020: Try Csaba Pinter's method for displaying image selected from dropdown menu or checkbox
        #I tried it and it does in fact allow you to switch images while retaining the zoom and slice position chosen
        #by the user in the prior image
        sliceCompositeNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
        for sliceCompositeNode in sliceCompositeNodes:
          sliceCompositeNode.SetBackgroundVolumeID(self.inputToMIP.GetID())
  
        #Return to side by side view that shows axial and sagittal views
        slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)

      slicer.mrmlScene.RemoveNode(self.sag_mip_node)

      self.inputSelector.setCurrentNode(self.inputToMIP)
      self.goToROICenter

      #Enable slice browsing with mouse and display slider for yellow slice
      #Enable mouse browsing
      #interactorStyle = slicer.app.layoutManager().sliceWidget('Yellow').sliceView().sliceViewInteractorStyle()
      #interactorStyle.SetActionEnabled(interactorStyle.BrowseSlice, True)
      #Add red slice slider to view
      lm = slicer.app.layoutManager()
      lm.sliceWidget('Yellow').sliceController().setVisible(True)
      
      #slicer.util.resetSliceViews() #want to always reset slider positions to center of image after going from MIP to non-MIP

      #7/21/2021: re-enabling this to see what happens
##      self.setAxialSliderPositionToROICenter(self.inputSelector.currentNode(),self.sag_mip_disp) #set sagittal position to center of ROI




  #Function for creating and showing axial MIP using method shared by Andrey
##  def showAxialMIP(self):
##    #first, get # of slices in image so you can use this for creating MIP
##    inputVolume = self.inputSelector.currentNode()
##    img = inputVolume.GetImageData()
##    rows,cols,slices = img.GetDimensions()
##
##    #Setup for altering image in red slice
##    sliceNodeRed = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed')
##    appLogic = slicer.app.applicationLogic()
##    sliceLogicRed = appLogic.GetSliceLogic(sliceNodeRed)
##    sliceLayerLogicRed = sliceLogicRed.GetBackgroundLayer()
##    resliceRed = sliceLayerLogicRed.GetReslice()
##
##    #If axial MIP box is checked, compute and display axial MIP of image
##    if(self.axchkbox.isChecked() == True):
##      print("Showing axial MIP")
##
##      #create MIP with SlabNumberOfSlices equal to number of z slices in original image
##      resliceRed.SetSlabModeToMax()
##      resliceRed.SetSlabNumberOfSlices(slices) #use very large value to ensure only one MIP (for now)
##      sliceNodeRed.Modified()
##
##      #Use red slice only view
##      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)
##
##    #If axial MIP box is unchecked, revert back to original image in red slice
##    if(self.axchkbox.isChecked() == False):
##      print("Showing original axial image")
##      #revert to original image in red slice by setting SlabNumberOfSlices to 1
##      resliceRed.SetSlabNumberOfSlices(1)
##      sliceNodeRed.Modified()
##
##      #Return to side by side view that shows axial and sagittal views
##      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)
##
##
##      
##  #Function for creating and showing sagittal MIP using method shared by Andrey
##  def showSagittalMIP(self):
##    #first, get # of slices in image so you can use this for creating MIP
##    inputVolume = self.inputSelector.currentNode()
##    img = inputVolume.GetImageData()
##    rows,cols,slices = img.GetDimensions()
##
##    #Setup for altering image in yellow slice
##    sliceNodeYellow = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeYellow')
##    appLogic = slicer.app.applicationLogic()
##    sliceLogicYellow = appLogic.GetSliceLogic(sliceNodeYellow)
##    sliceLayerLogicYellow = sliceLogicYellow.GetBackgroundLayer()
##    resliceYellow = sliceLayerLogicYellow.GetReslice()
##
##    #If sagittal MIP is checked, compute and display sagittal MIP of image
##    if(self.sagchkbox.isChecked() == True):
##      print("Showing sagittal MIP")
##
##      #create MIP with SlabNumberOfSlices equal to half the number of x slices in original image
##      #figure out if this is correct way to have a left breast sagittal and a right breast sagittal
##      resliceYellow.SetSlabModeToMax()
##      resliceYellow.SetSlabNumberOfSlices(1200) #use very large value to ensure only one MIP (for now)
##      sliceNodeYellow.Modified()
##
##      #Use yellow slice view
##      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView)
##
##    #If sagittal MIP box is unchecked, revert back to original image in yellow slice
##    if(self.sagchkbox.isChecked() == False):
##      print("Showing original sagittal image")
##      #revert to original image in yellow slice by setting SlabNumberOfSlices to 1
##      resliceYellow.SetSlabNumberOfSlices(1)
##      sliceNodeYellow.Modified()
##
##      #Return to side by side view that shows axial and sagittal views
##      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)

  
  def onSelect(self):
    #Comment this out because it makes Slicer crash
    #Edit 7/9/2020: Code for removing subraction MIP nodes
##    if ('subtraction axial MIP' in self.currentnodename):
##      #if you checked subtraction box 1st then ax MIP box, remove ax MIP node. Otherwise, remove subimg_node
##      if('subtraction' in self.inputToMIP.GetName()):
##        slicer.mrmlScene.RemoveNode(self.ax_mip_node)
##      else:
##        slicer.mrmlScene.RemoveNode(self.subimg_node)
##        
##    if ('subtraction sagittal MIP' in self.currentnodename):
##      #if you checked subtraction box 1st then sag MIP box, remove sag MIP node. Otherwise, remove subimg_node
##      if('subtraction' in self.inputToMIP.GetName()):
##        slicer.mrmlScene.RemoveNode(self.sag_mip_node)
##      else:
##        slicer.mrmlScene.RemoveNode(self.subimg_node)

    self.switchimage = True #variable to prevent calling of adjust window level function while changing image
                            #displayed from dropdown menu

    #Edit 9/29/2020:
    #If Module 2 is just starting (there is no startInputVolume), take window & level from current image.
    
    #Otherwise, take values from prior image.
    #Get window-level values set by user with
    #Adjust Window Level tool and Set the same values for the
    #next image in the dropdown menu (or from checkbox) that
    #the user switches to

    if(hasattr(self,'startInputVolume') == 0):
      self.window = self.inputSelector.currentNode().GetDisplayNode().GetWindow()
      self.level = 0.5*(self.inputSelector.currentNode().GetDisplayNode().GetWindowLevelMin() + self.inputSelector.currentNode().GetDisplayNode().GetWindowLevelMax())
    else:
      try:
        self.window = self.startInputVolume.GetDisplayNode().GetWindow()
        self.level = 0.5*(self.startInputVolume.GetDisplayNode().GetWindowLevelMin() + self.startInputVolume.GetDisplayNode().GetWindowLevelMax())
      except:
        #7/20/2021: If error, input volume was a subtraction
        self.window = self.subwindow
        self.level = self.sublevel
        
    self.currentnodename = self.inputSelector.currentNode().GetName()
    #If not currently viewing MIP, uncheck box for axial & sagittal MIPs
    if('MIP' not in self.currentnodename):
      self.newMIPSelection = 1 #variable to prevent overriding of user selection with MIP input after MIP box is unchecked

    if('subtraction' not in self.currentnodename):
      self.newSubSelection = 1 #variable to prevent overriding of user selected pre-contrast image with subtraction input


      #self.axchkbox.setChecked(False)
      #self.sagchkbox.setChecked(False)

    #First, display whatever user selected from dropdown menu
    self.applyButton.enabled = self.inputSelector.currentNode()
    #5/14/2020: Confirmed that this is how you make image shown in Slicer change to match user selection from drop down menu
    #slicer.util.setSliceViewerLayers(background=self.inputSelector.currentNodeID)

    #Edit 10/1/2020: Try Csaba Pinter's method for displaying image selected from dropdown menu
    #I tried it and it does in fact allow you to switch images while retaining the zoom and slice position chosen
    #by the user in the prior image
    sliceCompositeNodes = slicer.util.getNodesByClass('vtkMRMLSliceCompositeNode')
    for sliceCompositeNode in sliceCompositeNodes:
      sliceCompositeNode.SetBackgroundVolumeID(self.inputSelector.currentNode().GetID())



    #Edit 7/10/2020: If it exists, make segment node the foreground
    #Edit 8/24/2020: Only do this if there is a segment node and the SER colorization checkbox is checked
    if(hasattr(self,'segment_node') and self.segmentLesionCheckBox.isChecked() == True):
      slicer.util.setSliceViewerLayers(foreground=self.segment_node,foregroundOpacity = self.fopac)



    inputVolume = self.inputSelector.currentNode()
    inputVolume.GetDisplayNode().SetAutoWindowLevel(False) #Do this to prevent auto window leveling
    inputVolume.GetDisplayNode().SetWindowLevel(self.window,self.level) #set window level to match prior image
    img = inputVolume.GetImageData()
    rows,cols,slices = img.GetDimensions()
    #read inputVolume to numpy array
    sc = img.GetPointData().GetScalars()
    img_np = vtk_to_numpy(sc)
    #I tried other dimension orders but they gave "Scrambled" image
    img_np = img_np.reshape(slices,rows,cols) #3/24/2020: apparently this is the correct way to reshape numpy array so that it can be viewed in slicer as NIFTI
    img_np = img_np.astype('float64')
##    self.setSagittalSliderPositionToROICenter(inputVolume,img_np)


    #Edit 7/2/2020: Disabling this part because we eliminated dropdown menu method of MIP selection in favor of checkbox method
    
    #if axial MIP, adjust Slicer to only show red slice
##    if('MIPax' in self.currentnodename):
##      slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)
##    else:
##      #if sagittal MIP, adjust Slicer to only show yellow slice
##      if('MIPsag' in self.currentnodename):
##        slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView)
##      #If not a MIP, use side by side view that only shows axial and sagittal views
##      else:
##        slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutSideBySideView)

    #slicer.util.resetSliceViews() #get rid of this because this negates slice position settings for MIP

    #Only enable the subtraction checkbox if user is not viewing a pre-contrast image
    if('pre' in self.currentnodename):
      if(self.subtractCheckBox.isChecked() == True):
        slicer.mrmlScene.RemoveNode(self.subimg_node)
        self.subtractCheckBox.setChecked(False)
        
      self.subtractCheckBox.setEnabled(False)
      
    else:
      self.subtractCheckBox.setEnabled(True)

      #if subtraction box is checked but current user selected image is not a subtraction,
      #compute the subtraction of user selection from dropdown menu and display it
      if(self.subtractCheckBox.isChecked() == True and 'subtraction' not in self.currentnodename and self.axchkbox.isChecked() == False and self.axchkbox.isChecked() == False):
        slicer.mrmlScene.RemoveNode(self.subimg_node) #first, remove node for prior subtraction image that you were viewing
        self.showSubtractionFromNode()
          

    #The following 2 segments of code for MIPs are analogous to what I did for subtractions in the 3 lines above
    #These should be run regardless of whether user selection is pre or post contrast
    
    #if axial MIP box is checked but current user selected image is not an axial MIP,
    #compute the axial MIP of user selection from dropdown menu and display it
    if(self.axchkbox.isChecked() == True and 'axial MIP' not in self.currentnodename):
      try:
        slicer.mrmlScene.RemoveNode(self.ax_mip_node) #first, remove node for prior axial MIP image that you were viewing
      except:
        print("nothing to remove")
        
      self.showAxialMIPFromNode()
      self.currentnodename = self.inputSelector.currentNode().GetName()

      #In addition, if subtraction box is checked, convert axial MIP to subtraction axial MIP
      if(self.subtractCheckBox.isChecked() == True):
        self.showSubtractionFromNode()
        slicer.mrmlScene.RemoveNode(self.ax_mip_node)

      

    #if sagittal MIP box is checked but current user selected image is not a sagittal MIP,
    #compute the sagittal MIP of user selection from dropdown menu and display it
    if(self.sagchkbox.isChecked() == True and 'sagittal MIP' not in self.currentnodename):
      try:
        slicer.mrmlScene.RemoveNode(self.sag_mip_node) #first, remove node for prior sagittal MIP image that you were viewing
      except:
        print("nothing to remove")
      
      self.showSagittalMIPFromNode()
      self.currentnodename = self.inputSelector.currentNode().GetName()

      #In addition, if subtraction box is checked, convert sagittal MIP to subtraction sagittal MIP
      if(self.subtractCheckBox.isChecked() == True):
        self.showSubtractionFromNode()
        slicer.mrmlScene.RemoveNode(self.sag_mip_node)

    #Edit 9/29/2020: David and Jessica do not want Auto Window Level every time displayed image is switched

    #Revert back to AutoWindowLevel and set values in the sliders accordingly
##    inputVolume.GetDisplayNode().SetAutoWindowLevel(True) #Revert back to auto window level whenever you switch images
##    self.slider1.value = ( float( inputVolume.GetDisplayNode().GetWindow() ) ) 
##    self.slider2.value = ( float( inputVolume.GetDisplayNode().GetLevel() ) )


    self.switchimage = False #now that done switching image, set this to False


    #Edit 9/4/2020: Add lines of code to update text on upper left of slice views every time you switch image
    viewax = slicer.app.layoutManager().sliceWidget('Red').sliceView()
    if(self.axchkbox.isChecked() == True):
      viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,self.axmipnodename)
    else:
      viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,self.inputSelector.currentNode().GetName())

    viewsag = slicer.app.layoutManager().sliceWidget('Yellow').sliceView()
    if(self.sagchkbox.isChecked() == True):
      viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,self.sagmipnodename)
    else:
      viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperLeft,self.inputSelector.currentNode().GetName())

    #Edit 9/30/2020: Set startInputVolume so that you can access its window level settings and apply them to the
    #next image selected from dropdown menu
    self.startInputVolume = inputVolume


    
##    self.dfltrange = float(self.disprange.value)
##    self.dfltctr = float(self.dispctr.value)

    #Edit 8/11/2020: test out syntax for setting window level
    #https://radiopaedia.org/articles/windowing-ct?lang=us
    #the above link suggests that window width is the display range and
    #window level is center value in that range.
    #Therefore, when I tested this out with window = 10 and level = 1231,
    #it displayed the pre-contrast image with min = 1226 and max = 1236
    #inputVolume.GetScalarVolumeDisplayNode().SetAutoWindowLevel(False)
    #inputVolume.GetScalarVolumeDisplayNode().SetWindowLevel(10,1231)

    





  def onApplyButton(self):
    #Edit 5/21/2020: take absolute value of radius coordinates to prevent
    #errors from negative values
    #Specifically, was getting this error -- IndexError: cannot do a non-empty take from an empty axes.


    inputVolume = self.inputSelector.currentNode()
    
    #first, add final coordinates of ROI to arrays that belong to widget class object
    #these are already in voxel coordinates, all you have to do is convert them to type int
    self.roi.GetXYZ(self.roicenter[0,:])
    
    self.roi.GetRadiusXYZ(self.roiradius[0,:])
    self.roiradius = np.absolute(self.roiradius)

    #convert ROI coordinates from RAS to IJK
    self.roicenter[0,:], self.roiradius[0,:] = RASToIJKFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)

    print("roi center")
    print(self.roicenter)
    print("roi radius")
    print(self.roiradius)

    #Edit 7/30/2020: Pressing "Generate Report" button will no longer trigger removal of vtkMRMLAnnotationROINode's
    #Remove ROI from scene
##    slicer.mrmlScene.RemoveNode(self.roi)

    #Tried the line below and got error like "self.roi has no attribute WriteXML"
    #self.roi.WriteXML("ROI.xml",1)

    #then, add final coordinates of all omit regions to arrays that belong to widget class object
    #these are already in voxel coordinates, all you have to do is convert them to type int
    self.omit1.GetXYZ(self.omitcenters[0,:])
    self.omit2.GetXYZ(self.omitcenters[1,:])
    self.omit3.GetXYZ(self.omitcenters[2,:])
    self.omit4.GetXYZ(self.omitcenters[3,:])
    self.omit5.GetXYZ(self.omitcenters[4,:])

    self.omit1.GetRadiusXYZ(self.omitradii[0,:])
    self.omit2.GetRadiusXYZ(self.omitradii[1,:])
    self.omit3.GetRadiusXYZ(self.omitradii[2,:])
    self.omit4.GetRadiusXYZ(self.omitradii[3,:])
    self.omit5.GetRadiusXYZ(self.omitradii[4,:])
    self.omitradii = np.absolute(self.omitradii)

    #convert all omit coordinates from RAS to IJK
    self.omitcenters[0,:], self.omitradii[0,:] = RASToIJKFunc(self.omitcenters[0,:],self.omitradii[0,:],inputVolume)
    self.omitcenters[1,:], self.omitradii[1,:] = RASToIJKFunc(self.omitcenters[1,:],self.omitradii[1,:],inputVolume)
    self.omitcenters[2,:], self.omitradii[2,:] = RASToIJKFunc(self.omitcenters[2,:],self.omitradii[2,:],inputVolume)
    self.omitcenters[3,:], self.omitradii[3,:] = RASToIJKFunc(self.omitcenters[3,:],self.omitradii[3,:],inputVolume)
    self.omitcenters[4,:], self.omitradii[4,:] = RASToIJKFunc(self.omitcenters[4,:],self.omitradii[4,:],inputVolume)

    
    print("omit centers")
    print(self.omitcenters)

    print("omit radii")
    print(self.omitradii)

    #Edit 7/30/2020: Pressing "Generate Report" button will no longer trigger removal of vtkMRMLAnnotationROINode's
    #Remove all omits defined by user from scene
##    for i in range(self.omitCount):
##      if i == 0:
##        slicer.mrmlScene.RemoveNode(self.omit1)
##        
##      if i == 1:
##        slicer.mrmlScene.RemoveNode(self.omit2)
##        
##      if i == 2:
##        slicer.mrmlScene.RemoveNode(self.omit3)
##        
##      if i == 3:
##        slicer.mrmlScene.RemoveNode(self.omit4)
##        
##      if i == 4:
##        slicer.mrmlScene.RemoveNode(self.omit5)
        
    #Edit 6/29/2020: Compare z value of intial position in IJKToRAS matrix from inputVolume to same element in
    #my IJKToLPS aff_mat, and if they do match, add an extra flip to sagittal images in report
    #If these two z origin values match, this means Slicer used closest to foot as 1st slice and
    #closest to head as last slice, which is why you have to flip
    ijkToRASmat = vtk.vtkMatrix4x4()
    inputVolume.GetIJKToRASMatrix(ijkToRASmat)

    #Edit 9/29/2020:
    #If Module 2 is just starting (there is no startInputVolume), take window & level from current image.
    
    #Otherwise, take values from prior image.

    #Either way, these values will be used for report
    
    if(hasattr(self,'startInputVolume') == 0):
      self.window = self.inputSelector.currentNode().GetDisplayNode().GetWindow()
      self.level = 0.5*(self.inputSelector.currentNode().GetDisplayNode().GetWindowLevelMin() + self.inputSelector.currentNode().GetDisplayNode().GetWindowLevelMax())
    else:
      self.window = self.startInputVolume.GetDisplayNode().GetWindow()
      self.level = 0.5*(self.startInputVolume.GetDisplayNode().GetWindowLevelMin() + self.startInputVolume.GetDisplayNode().GetWindowLevelMax())

    #Edit 9/30/2020: For report, window and level should always be adjusted upward compared to what is used for displaying image in Slicer
    #Still getting black images for some exams. Try multiplying by 2 and 2.5 instead of adding 2500 and 3000
    window_report = self.window*2
    level_report = self.window*2.5

    #Edit 10/6/2020: retrieve threshold values from spin boxes and use them for tumor mask calculations for report
    print("threshold values")
    #NOTE: The way this code is set up, the BKG thresh value from spin box MUST BE DIVIDED BY 100
    pct = self.bkgthresh.value/100 
    print(pct)
    pethresh = self.pethresh.value
    print(pethresh)
    minconnpix = self.minconnthresh.value
    print(minconnpix)

    #Edit 6/30/2020: no more imgReversed
    #Start the "run" section processing, using the final ROI selections
    logic = DCE_TumorMapProcessLogic()
    logic.run(self.gzipped,self.exampath,self.inputSelector.currentNode(),self.roiradius,self.roicenter,self.omitCount,self.omitradii,self.omitcenters,ijkToRASmat,self.tempres, self.fsort, self.nslice, self.manufacturer, self.studydate, self.dce_folders, self.aff_mat, self.aff_inv_mat, self.earlyPostContrastNum, self.latePostContrastNum, self.earlydiffmm, self.earlydiffss, self.latediffmm, self.latediffss,self.sitestr,self.idstr,self.nodevisstr,window_report,level_report,pct,pethresh,minconnpix)

  def addOmitRegion(self):
    
    self.omitCount = self.omitCount+1

    #Edit 6/8/2020: Omit always initialized at center of image regardless of image size
    inputVolume = self.inputSelector.currentNode()
    img = inputVolume.GetImageData()
    rows,cols,slices = img.GetDimensions()

    roicenter = np.array([int(rows/2),int(cols/2),int(slices/2)])
    #Edit 7/19/2021: Give omit region 25% smaller initial size to avoid
    #confusion with the ROI.
    roiradius = np.array([100,100,40])

    #Call function that returns RAS coordinate center and radius give input volume and IJK coordinate center and radius
    roicenter_RAS, roiradius_RAS = IJKToRASFunc(roicenter,roiradius,inputVolume)
    #7/19/2021: Initial omit radius should be 75% of initial roiradius
    omitradius_RAS = 0.75*roiradius_RAS
    #7/16/2021:Add yes/no dialog to this to avoid accidental omit region additions.
    if(self.omitCount <= 5):
      yesnotxt = "Do you want to add omit region #" + str(self.omitCount) + "? If yes, omit region will be placed at center of image."
    else:
      slicer.util.confirmOkCancelDisplay("You have already added the maximum allowable number of omit regions.","Cannot add another omit region.")

      
    #Choose which omit object to add to scene based on omit count
    if self.omitCount == 1:
      if(slicer.util.confirmYesNoDisplay(yesnotxt)):
        self.omit1.SetXYZ(roicenter_RAS)
        self.omit1.SetRadiusXYZ(omitradius_RAS)
        slicer.mrmlScene.AddNode(self.omit1)
        slicer.util.resetSliceViews()
      else:
        self.omitCount = self.omitCount - 1
        
    if self.omitCount == 2:
      if(slicer.util.confirmYesNoDisplay(yesnotxt)):
        self.omit2.SetXYZ(roicenter_RAS)
        self.omit2.SetRadiusXYZ(omitradius_RAS)
        slicer.mrmlScene.AddNode(self.omit2)
        slicer.util.resetSliceViews()
      else:
        self.omitCount = self.omitCount - 1

      
    if self.omitCount == 3:
      if(slicer.util.confirmYesNoDisplay(yesnotxt)):
        self.omit3.SetXYZ(roicenter_RAS)
        self.omit3.SetRadiusXYZ(omitradius_RAS)
        slicer.mrmlScene.AddNode(self.omit3)
        slicer.util.resetSliceViews()
      else:
        self.omitCount = self.omitCount - 1

      
    if self.omitCount == 4:
      if(slicer.util.confirmYesNoDisplay(yesnotxt)):
        self.omit4.SetXYZ(roicenter_RAS)
        self.omit4.SetRadiusXYZ(omitradius_RAS)
        slicer.mrmlScene.AddNode(self.omit4)
        slicer.util.resetSliceViews()
      else:
        self.omitCount = self.omitCount - 1

      
    if self.omitCount == 5:
      if(slicer.util.confirmYesNoDisplay(yesnotxt)):
        self.omit5.SetXYZ(roicenter_RAS)
        self.omit5.SetRadiusXYZ(omitradius_RAS)
        slicer.mrmlScene.AddNode(self.omit5)
        slicer.util.resetSliceViews()
      else:
        self.omitCount = self.omitCount - 1

      
    if self.omitCount > 5:
      print("Cannot add any more omit regions")


  def importROIFromFile(self):
    #Edit 6/26/2020: Instead of doing LPS --> IJK --> RAS and using that to set ROI position,
    #I will see if I can just take the LPS coordinates, multiply the first two values by -1,
    #and input that directly into node position setting function. I'm trying this because both
    #David and Andrey said that slice order should not affect the ability of the RAS (or LPS)
    #coordinates in xml file to correspond to tumor
    
##    pre_folder = int(dce_folders[0])
##    print("Pre-contrast folder # is: " + str(pre_folder))


    #Edit 9/14/2020: message box interface for choosing ROI xml file for any visit for current patient

    if('ispy' in self.exampath or 'acrin' in self.exampath):
      #find all visit folders to come up with button labels for message box
      folders = os.listdir(self.idpath)
      #Edit 2/4/2021: restrict folders to only include folders with 'v' in the name
      folders = [f for f in folders if 'v' in f]
      print(folders)
      visit_strs = []
      for f in folders:
        vpos = f.find('v')
        vis_id = int(f[vpos+1:vpos+3])

        #date is either the 6 digit number just before the '_' or just after the '_'
        try:
          vis_date = f.split('_')[0]
          
          #Edit 12/14/2020: if statement for case when this split gives you visit instead of date
          if('v' in vis_date):
            vis_date = f.split('_')[1]
            
        except:
          vis_date = f.split('_')[1]


        if(vis_id%10 == 0):
          visstr = 'MR' + str(int(vis_id/10)) + '\n' + vis_date
        else:
          visstr = 'MR' + str(float(vis_id/10)) + '\n' + vis_date

        visit_strs.append(visstr)

      print(visit_strs)

      if len(visit_strs) > 1:
        #Message box that allows user to choose which visit to load ROI for?
        choosevisbox = qt.QMessageBox()
        choosevisbox.setStandardButtons(0) #remove default "Ok" button
        choosevisbox.setText('Which visit''s ROI do you want to load?')

        #loop for adding buttons corresponding to visits
        for vs in range(len(visit_strs)):
          if(vs == 0):
            button1 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)
          
          if(vs == 1):
            button2 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)
          
          if(vs == 2):
            button3 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)
          
          if(vs == 3):
            button4 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)
          
          if(vs == 4):
            button5 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)
          
          if(vs == 5):
            button6 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)
          
          if(vs == 6):
            button7 = choosevisbox.addButton(visit_strs[vs],qt.QMessageBox.YesRole)

        #show the choose visit box
        choosevisbox.exec_()

        #set the default path for qfiledialog xml file selection based on user selection from choose visit message box
        #doing this in a loop because it is not certain for me how many buttons there will be; it depends on which folder user selects
        for sel in range(len(folders)):
          if sel == 0:
            if(choosevisbox.clickedButton() == button1):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"

          if sel == 1:
            if(choosevisbox.clickedButton() == button2):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"

          if sel == 2:
            if(choosevisbox.clickedButton() == button3):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"
             
          if sel == 3:
            if(choosevisbox.clickedButton() == button4):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"
             
          if sel == 4:
            if(choosevisbox.clickedButton() == button5):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"
             
          if sel == 5:
            if(choosevisbox.clickedButton() == button6):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"
             
          if sel == 6:
            if(choosevisbox.clickedButton() == button7):
              visitpath = self.idpath + "\\" + folders[sel]
              examfolders = os.listdir(visitpath)
              exampath = visitpath + "\\" + examfolders[0]
              self.dflt_xml_dir = exampath + "\\" + "voi_lps_files"
              
      #If only one visit folder, choose that as your xml default directory
      else:
        self.dflt_xml_dir = self.exampath + "\\" + "voi_lps_files"
    else:
      #If not \\researchfiles MR exam directory structure
      #Edit 6/29/2021: David wants this to default to the
      #patient folder instead. Try that, and default to
      #exampath if you get an error.
      #find indices where slashes '/' occur in exampath
      slashinds = []
      for i in range(len(self.exampath)):
        if(self.exampath[i] == '/'):
          slashinds.append(i)
      try:
        self.dflt_xml_dir = self.exampath[0:slashinds[-2]+1]
      except:
        self.dflt_xml_dir = self.exampath

    #self.dflt_xml_dir = self.exampath + "\\" + "voi_lps_files"
    
    #allow user to select xml file they want to import ROI from, if the exam has existing xml files
    try:
      xmlfilepath = qt.QFileDialog.getOpenFileName(0,"Select xml file to import ROI from",self.dflt_xml_dir)
    except:
      print("no existing xml file to import")
      return

    roicenter, roiradius, ocenters, oradii = compute_lps_to_rcs.readRASCoordsFromXMLFile(xmlfilepath)

    #call function that returns voxel bounds of ROI from file
##    xs, xf, ys, yf, zs, zf, oxs, oxf, oys, oyf, ozs, ozf, voi_mask, aff_mat, aff_inv_mat = compute_lps_to_rcs.getVOIVoxelsFromInverseAffine(self.exampath,xmlfilepath,pre_folder,nslice,fsort)
##
##    inputVolume = self.inputSelector.currentNode()
##
##    roicenter = np.zeros((1,3))
##    roicenter[0,0] = int((xs+xf)/2)
##    roicenter[0,1] = int((ys+yf)/2)
##    roicenter[0,2] = int((zs+zf)/2)
##
##    roiradius = np.zeros((1,3))
##    roiradius[0,0] = int(roicenter[0,0])-xs
##    roiradius[0,1] = int(roicenter[0,1])-ys
##    roiradius[0,2] = int(roicenter[0,2])-zs
##    roiradius = np.abs(roiradius)
##
##    #Convert ROI center and radius to RAS before using to set AnnotationROINode position
##    roicenter, roiradius = IJKToRASFunc(roicenter,roiradius,inputVolume)
##
##    print(roicenter)
##    print(roiradius)


    #Remove ROI from scene and then add it back to scene with
    #position from xml file
    #slicer.mrmlScene.RemoveNode(self.roi)
    self.roi.SetXYZ(roicenter)
    self.roi.SetRadiusXYZ(roiradius)
    #slicer.mrmlScene.AddNode(self.roi)

    #7/16/2021: Before proceeding to check xml file for omits,
    #remove any omits that may already be in the scene.
    if(self.omitCount > 0):
      for oc in range(self.omitCount):
        if(oc == 0):
          slicer.mrmlScene.RemoveNode(self.omit1)

        if(oc == 1):
          slicer.mrmlScene.RemoveNode(self.omit2)

        if(oc == 2):
          slicer.mrmlScene.RemoveNode(self.omit3)

        if(oc == 3):
          slicer.mrmlScene.RemoveNode(self.omit4)

        if(oc == 4):
          slicer.mrmlScene.RemoveNode(self.omit5)
      self.omitCount = 0

    if ocenters[0,0] != -1:
      omits = 1
    else:
      omits = 0

    #Edit 6/26/2020: Add omits to ROI nodes using outputs from new RAS function
    if omits == 1:
      for i in range(ocenters.shape[0]):
        self.omitCount = self.omitCount+1

        if self.omitCount > 5:
          #Edit 2/2/21: Set omitCount to 5 here to avoid errors later on in code.
          self.omitCount = 5
          print("Cannot add any more omit regions")
          break

        curr_center = ocenters[i,:]
        curr_radius = oradii[i,:]
        
        #Choose which omit object to add to scene based on omit count
        if self.omitCount == 1:
          self.omit1.SetXYZ(curr_center)
          self.omit1.SetRadiusXYZ(curr_radius)
          slicer.mrmlScene.AddNode(self.omit1)
      
        if self.omitCount == 2:
          self.omit2.SetXYZ(curr_center)
          self.omit2.SetRadiusXYZ(curr_radius)
          slicer.mrmlScene.AddNode(self.omit2)
      
        if self.omitCount == 3:
          self.omit3.SetXYZ(curr_center)
          self.omit3.SetRadiusXYZ(curr_radius)
          slicer.mrmlScene.AddNode(self.omit3)
      
        if self.omitCount == 4:
          self.omit4.SetXYZ(curr_center)
          self.omit4.SetRadiusXYZ(curr_radius)
          slicer.mrmlScene.AddNode(self.omit4)
      
        if self.omitCount == 5:
          self.omit5.SetXYZ(curr_center)
          self.omit5.SetRadiusXYZ(curr_radius)
          slicer.mrmlScene.AddNode(self.omit5)

    #Edit 7/31/2020: Set red and yellow slice positions to center axial and sagittal slices of ROI
    inputVolume = self.inputSelector.currentNode()
    img = getNPImgFromNode(self.inputSelector.currentNode().GetName())
    slicer.util.resetSliceViews()
    self.setAxialSliderPositionToROICenter(inputVolume,img)
    self.setSagittalSliderPositionToROICenter(inputVolume,img)


    

    #Edit 5/21/2020: Import omits from file as well
    #Check if there are any omits at all
##    if oxs[0] != -1:
##      for i in range(len(oxs)):
##        self.omitCount = self.omitCount+1
##
##        if self.omitCount > 5:
##          print("Cannot add any more omit regions")
##          break
##
##        #create np arrays containing center and radius of ith omit region
##        curr_center = np.zeros((1,3))
##        curr_center[0,0] = ( int(oxs[i]) + int(oxf[i]) )/2
##        curr_center[0,1] = ( int(oys[i]) + int(oyf[i]) )/2
##        curr_center[0,2] = ( int(ozs[i]) + int(ozf[i]) )/2
##
##
##        curr_radius = np.zeros((1,3))
##        curr_radius[0,0] = abs( int(curr_center[0,0]) - int(oxs[i]) )
##        curr_radius[0,1] = abs( int(curr_center[0,1]) - int(oys[i]) )
##        curr_radius[0,2] = abs( int(curr_center[0,2]) - int(ozs[i]) )
##
##        
##        #Convert current omit's center and radius to RAS before using to set AnnotationROINode position
##        curr_center, curr_radius = IJKToRASFunc(curr_center,curr_radius,inputVolume)
##
##        #Choose which omit object to add to scene based on omit count
##        if self.omitCount == 1:
##          self.omit1.SetXYZ(curr_center)
##          self.omit1.SetRadiusXYZ(curr_radius)
##          slicer.mrmlScene.AddNode(self.omit1)
##      
##        if self.omitCount == 2:
##          self.omit2.SetXYZ(curr_center)
##          self.omit2.SetRadiusXYZ(curr_radius)
##          slicer.mrmlScene.AddNode(self.omit2)
##      
##        if self.omitCount == 3:
##          self.omit3.SetXYZ(curr_center)
##          self.omit3.SetRadiusXYZ(curr_radius)
##          slicer.mrmlScene.AddNode(self.omit3)
##      
##        if self.omitCount == 4:
##          self.omit4.SetXYZ(curr_center)
##          self.omit4.SetRadiusXYZ(curr_radius)
##          slicer.mrmlScene.AddNode(self.omit4)
##      
##        if self.omitCount == 5:
##          self.omit5.SetXYZ(curr_center)
##          self.omit5.SetRadiusXYZ(curr_radius)
##          slicer.mrmlScene.AddNode(self.omit5)



#Edit 6/30/2020: invertROIPos will be commented out because it is no longer needed
  #Edit 6/25/2020: Function to use in case where ROI doesn't match with tumor because image appears to have slices in reverse order
##  def invertROIpos(self):
##    self.imgReversed = 1 #This is how to track if image has axial slices reversed for report
##
##    inputVolume = self.inputSelector.currentNode()
##    img = inputVolume.GetImageData()
##    rows,cols,slices = img.GetDimensions()
##
##    #read ROI and omit region info
##    self.roi.GetXYZ(self.roicenter[0,:])
##    
##    self.roi.GetRadiusXYZ(self.roiradius[0,:])
##
##    #convert ROI coordinates from RAS to IJK
##    self.roicenter[0,:], self.roiradius[0,:] = RASToIJKFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)
##    #invert z position
##    self.roicenter[0,2] = slices - self.roicenter[0,2] 
##    #convert from IJK to RAS again
##    self.roicenter[0,:], self.roiradius[0,:] = IJKToRASFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)
##    #Add new z-position to node
##    self.roi.SetXYZ(self.roicenter[0,:])
##    self.roi.SetRadiusXYZ(self.roiradius[0,:])
##
##    #repeat the process above for all omits that are in the scene
##    for o in range(self.omitCount):
##      if o == 0:
##        self.omit1.GetXYZ(self.omitcenters[0,:])
##        self.omit1.GetRadiusXYZ(self.omitradii[0,:])
##        self.omitcenters[0,:], self.omitradii[0,:] = RASToIJKFunc(self.omitcenters[0,:],self.omitradii[0,:],inputVolume)
##        self.omitcenters[0,2] = slices - self.omitcenters[0,2]
##        self.omitcenters[0,:], self.omitradii[0,:] = IJKToRASFunc(self.omitcenters[0,:],self.omitradii[0,:],inputVolume)
##        self.omit1.SetXYZ(self.omitcenters[0,:])
##        self.omit1.SetRadiusXYZ(self.omitradii[0,:])
##
##      if o == 1:
##        self.omit2.GetXYZ(self.omitcenters[1,:])
##        self.omit2.GetRadiusXYZ(self.omitradii[1,:])
##        self.omitcenters[1,:], self.omitradii[1,:] = RASToIJKFunc(self.omitcenters[1,:],self.omitradii[1,:],inputVolume)
##        self.omitcenters[1,2] = slices - self.omitcenters[1,2]
##        self.omitcenters[1,:], self.omitradii[1,:] = IJKToRASFunc(self.omitcenters[1,:],self.omitradii[1,:],inputVolume)
##        self.omit2.SetXYZ(self.omitcenters[1,:])
##        self.omit2.SetRadiusXYZ(self.omitradii[1,:])
##
##      if o == 2:
##        self.omit3.GetXYZ(self.omitcenters[2,:])
##        self.omit3.GetRadiusXYZ(self.omitradii[2,:])
##        self.omitcenters[2,:], self.omitradii[2,:] = RASToIJKFunc(self.omitcenters[2,:],self.omitradii[2,:],inputVolume)
##        self.omitcenters[2,2] = slices - self.omitcenters[2,2]
##        self.omitcenters[2,:], self.omitradii[2,:] = IJKToRASFunc(self.omitcenters[2,:],self.omitradii[2,:],inputVolume)
##        self.omit3.SetXYZ(self.omitcenters[2,:])
##        self.omit3.SetRadiusXYZ(self.omitradii[2,:])
##
##      if o == 3:
##        self.omit4.GetXYZ(self.omitcenters[3,:])
##        self.omit4.GetRadiusXYZ(self.omitradii[3,:])
##        self.omitcenters[3,:], self.omitradii[3,:] = RASToIJKFunc(self.omitcenters[3,:],self.omitradii[3,:],inputVolume)
##        self.omitcenters[3,2] = slices - self.omitcenters[3,2]
##        self.omitcenters[3,:], self.omitradii[3,:] = IJKToRASFunc(self.omitcenters[3,:],self.omitradii[3,:],inputVolume)
##        self.omit4.SetXYZ(self.omitcenters[3,:])
##        self.omit4.SetRadiusXYZ(self.omitradii[3,:])
##
##      if o == 4:
##        self.omit5.GetXYZ(self.omitcenters[4,:])
##        self.omit5.GetRadiusXYZ(self.omitradii[4,:])
##        self.omitcenters[4,:], self.omitradii[4,:] = RASToIJKFunc(self.omitcenters[4,:],self.omitradii[4,:],inputVolume)
##        self.omitcenters[4,2] = slices - self.omitcenters[4,2]
##        self.omitcenters[4,:], self.omitradii[4,:] = IJKToRASFunc(self.omitcenters[4,:],self.omitradii[4,:],inputVolume)
##        self.omit5.SetXYZ(self.omitcenters[4,:])
##        self.omit5.SetRadiusXYZ(self.omitradii[4,:])



  def segmentLesion(self):
    #Edit 5/21/2020: take absolute value of radius coordinates to prevent
    #errors from negative values
    #Specifically, was getting this error -- IndexError: cannot do a non-empty take from an empty axes.

    self.fopac = 0.4 #Edit 10/1/2020: use this variable for foreground opacity

    #Edit 7/24/2020: If segmentLesion box is checked, compute and display SER colorized lesion segmentation
    if(self.segmentLesionCheckBox.isChecked() == True):

      #read ROI and omit region info
      self.roi.GetXYZ(self.roicenter[0,:])
    
      self.roi.GetRadiusXYZ(self.roiradius[0,:])
      self.roiradius = np.absolute(self.roiradius)

      print(self.roiradius)

      self.omit1.GetXYZ(self.omitcenters[0,:])
      self.omit2.GetXYZ(self.omitcenters[1,:])
      self.omit3.GetXYZ(self.omitcenters[2,:])
      self.omit4.GetXYZ(self.omitcenters[3,:])
      self.omit5.GetXYZ(self.omitcenters[4,:])

      self.omit1.GetRadiusXYZ(self.omitradii[0,:])
      self.omit2.GetRadiusXYZ(self.omitradii[1,:])
      self.omit3.GetRadiusXYZ(self.omitradii[2,:])
      self.omit4.GetRadiusXYZ(self.omitradii[3,:])
      self.omit5.GetRadiusXYZ(self.omitradii[4,:])
      self.omitradii = np.absolute(self.omitradii)

      inputVolume = self.inputSelector.currentNode()
      self.startInputVolume = inputVolume #volume to use for window leveling
      #Edit 8/21/2020: code for evaluating if you should use existing segment node or make new one
      use_exist = 0
      if(hasattr(self,'segment_node')):

        roicentercp = np.round(self.roicenter)
        segroicentercp = np.round(self.segroicenter)
        roicentermatch = (roicentercp == segroicentercp)
        print(roicentermatch)

        roiradiuscp = np.round(self.roiradius)
        print(roiradiuscp)
        segroiradiuscp = np.round(self.segroiradius)
        print(segroiradiuscp)
        roiradiusmatch = (roiradiuscp == segroiradiuscp)
        print(roiradiusmatch)
        print(roiradiuscp-segroiradiuscp)

        omitcenterscp = np.round(self.omitcenters)
        segomitcenterscp = np.round(self.segomitcenters)
        omitcentersmatch = (omitcenterscp == segomitcenterscp)
        print(omitcentersmatch)

        omitradiicp = np.round(self.omitradii)
        segomitradiicp = np.round(self.segomitradii)
        omitradiimatch = (omitradiicp == segomitradiicp)
        print(omitradiimatch)

        
        #only use existing if there is a segment node and current ROI and omit positions are same as that one
        #Edit 10/6/2020: ROI and all 3 thresholds (BKG, PE, MC) must be unchanged in order to reuse existing segmentation
        if(roicentermatch.all() and roiradiusmatch.all() and omitradiimatch.all() and omitcentersmatch.all() and self.bkgthresh.value/100 == self.segbkgthresh and self.pethresh.value == self.segpethresh and self.minconnthresh.value == self.segmcthresh):
          use_exist = 1
        

      #Edit 8/20/2020
      #if there is an existing segment node and the ROI and omit coordinates have not changed, load the existing segment node 
      if(use_exist == 1):
        print("using existing SER colorization")
        #slicer.mrmlScene.RemoveNode(inputVolume)
        
        slicer.util.setSliceViewerLayers(foreground=self.segment_node,foregroundOpacity = self.fopac)

        #slicer.util.setSliceViewerLayers(background=inputVolume)
        #slicer.util.resetSliceViews()

      #Otherwise, compute new SER colorization and display that
      else:
        print("making new SER colorization")

        #Add progress bar to update user on which step in process module is at
        progressBar = slicer.util.createProgressDialog(windowTitle = "New SER Colorization")
      
        img = inputVolume.GetImageData()
        rows,cols,slices = img.GetDimensions()
        print("rows: " + str(rows) + " cols: " + str(cols) + " slices: " + str(slices))

        #read inputVolume to numpy array
        sc = img.GetPointData().GetScalars()
        img_np = vtk_to_numpy(sc)
        #I tried other dimension orders but they gave "Scrambled" image

        img_np = img_np.reshape(slices,rows,cols) #3/24/2020: apparently this is the correct way to reshape numpy array so that it can be viewed in slicer as NIFTI
    ##    img_np = img_np.astype('float64') # convert to 64-bit datatype to avoid too small dynamic range 
    ##
        #normalization before making rgb appears to work well

        img_np = img_np.transpose(2,1,0) #transpose to cols,rows,slices to give same orientation as input. Can't do this earlier because image gets "scrambled"
        #7/13/2021: Can make int8 to avoid memory error
        #because only using this to retrieve shape.
        img_np = img_np.astype('int8')
        #img_np = 128*img_np/(np.amax(img_np)) #normalize to max 128
        #rgb version of input image that we will add SER colorization to

        #normalization sets everything to 0 no matter how I do it
    ##    img_np_rgb = img_np_rgb.astype('float64') #do this so dividing by max doesn't make everything 0
    ##    norm_array = np.amax(img_np_rgb)*np.ones(img_np_rgb.shape)
    ##    norm_array = norm_array.astype('float64')
    ##    img_np_rgb = img_np_rgb/norm_array

        #Update progress bar to say that SER colorized image array has been initialized
        progressBar.value = 20
        progressBar.labelText = "Initialized array for SER colorized image"
        slicer.app.processEvents()


        #since negative value voxels end up really bright in rgb image, set these to 0 first
##        img_np = np.where(img_np>0,img_np,0) #set to 0 where img_np is NOT greater than 0
##        img_np_u1 = 255*img_np/(np.amax(img_np)) #convert to 8-bit unsigned (u1) then add to rgb array



        #Edit 7/31/2020: img_np_rgb should be initialized as all zeros, NOT as input image repeated across all 3 color channels

        #define contents of array after datatype conversion
  ##      img_np_rgb[:,:,:,0] = img_np_u1
  ##      img_np_rgb[:,:,:,1] = img_np_u1
  ##      img_np_rgb[:,:,:,2] = img_np_u1

        #SER colorization -- need to read pre, early, and late images into numpy arrays, then compute PE and SER images, then make tumor mask
        #then use tumor mask and SER image to define colorized overlay

        #Using IJKToRAS matrix from inputVolume for nii is giving syntax errors, so try obtaining RAS matrix from LPS aff_mat
    
        #Save RAS coordinates of ROI and omits to new variables
        #These can be used to check if new lesion segmentations is needed or
        #if you can use existing one
        self.segroicenter = self.roicenter
        self.segroiradius = np.absolute(self.roiradius)

        self.segomitcenters = self.omitcenters
        self.segomitradii = np.absolute(self.omitradii)

        #convert ROI coordinates from RAS to IJK
        self.roicenter[0,:], self.roiradius[0,:] = RASToIJKFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)

        #convert all omit coordinates from RAS to IJK
        self.omitcenters[0,:], self.omitradii[0,:] = RASToIJKFunc(self.omitcenters[0,:],self.omitradii[0,:],inputVolume)
        self.omitcenters[1,:], self.omitradii[1,:] = RASToIJKFunc(self.omitcenters[1,:],self.omitradii[1,:],inputVolume)
        self.omitcenters[2,:], self.omitradii[2,:] = RASToIJKFunc(self.omitcenters[2,:],self.omitradii[2,:],inputVolume)
        self.omitcenters[3,:], self.omitradii[3,:] = RASToIJKFunc(self.omitcenters[3,:],self.omitradii[3,:],inputVolume)
        self.omitcenters[4,:], self.omitradii[4,:] = RASToIJKFunc(self.omitcenters[4,:],self.omitradii[4,:],inputVolume)

        #Update progress bar to say that reading of ROI and omit coordinates is done
        progressBar.value = 40
        progressBar.labelText = "Done reading ROI and omit voxel coordinates"
        slicer.app.processEvents()


        #generate FTV processing outputs
        #NOTE: The way this code is set up, the BKG thresh value from spin box MUST BE DIVIDED BY 100

        #threshold values used for lesion segmentation
        self.segbkgthresh = self.bkgthresh.value/100
        self.segpethresh = self.pethresh.value
        self.segmcthresh = self.minconnthresh.value

        pre_folder = int(self.dce_folders[0])
        #Edit 7/20/2020: If exam was originally gzipped, use imagespath as input to this function instead of exampath
        if(self.gzipped == 1):
          a,b,c,pe,self.ser,tumor_mask,voi_mask,zs,zf,ys,yf,xs,xf,pct,pre_thresh,pethresh,minconnpix = ftv_map_gen.makeFTVMaps(self.imagespath, self.manufacturer, self.dce_folders,self.roicenter,self.roiradius,self.omitcenters,self.omitradii, self.earlyPostContrastNum, self.latePostContrastNum,self.segbkgthresh,self.segpethresh,self.segmcthresh)
          #DICOM header path for pre-contrast image
          pre_path = self.imagespath + "\\" + str(pre_folder)
        #Otherwise, use exampath
        else:
          a,b,c,pe,self.ser,tumor_mask,voi_mask,zs,zf,ys,yf,xs,xf,pct,pre_thresh,pethresh,minconnpix = ftv_map_gen.makeFTVMaps(self.exampath, self.manufacturer, self.dce_folders,self.roicenter,self.roiradius,self.omitcenters,self.omitradii, self.earlyPostContrastNum, self.latePostContrastNum,self.segbkgthresh,self.segpethresh,self.segmcthresh)
          #DICOM header path for pre-contrast image
          pre_path = self.exampath + "\\" + str(pre_folder)


        #Update progress bar to say that SER map has been generated
        progressBar.value = 60
        progressBar.labelText = "SER map created"
        slicer.app.processEvents()


        #Edit 5/21/2020: Want to print ftv in cc every time user makes lesion segmentation
        #First, compute ftv in voxels
        #7/13/2021: Only multiply by voi_mask if maps are not cropped.
##        if(voi_mask.shape[0] == a.shape[0]):
        tumor_mask_voi = tumor_mask*voi_mask
##        else:
##          tumor_mask_voi = tumor_mask
          
        ftv_voxels = np.sum(tumor_mask_voi) #ftv in # of voxels
        print("Number of voxels in tumor mask")
        print(ftv_voxels)

        #4/23/21: Print number of voxels with SER 0 to see if not applying this
        #filter is significantly affecting results
##        if(self.ser.shape[0] == tumor_mask_voi.shape[0]):
##          ser_masked = self.ser*tumor_mask_voi
##        else:
##          serfull = np.zeros((tumor_mask_voi.shape))
##          serfull[xs:xf+1,ys:yf+1,zs:zf+1] = self.ser
##          ser_masked = self.ser*tumor_mask_voi
##          
##        ser_maskednegval = np.sum(ser_masked < 0)
##        print("Number of VOI voxels with SER less than 0")
##        print(ser_maskednegval)

    
        #Then, use DICOM header to compute voxel volume in cubic centimeters (cc)
        pre_imgs = os.listdir(pre_path)
        pre_imgs = sorted(pre_imgs)
        pre_img1path = pre_path + "\\" + pre_imgs[0]
    
        try:
          pre_hdr1 = pydicom.dcmread(pre_img1path,stop_before_pixels = True)
        except:
          pre_hdr1 = dicom.read_file(pre_img1path)

        voxsize_mm3 = float(pre_hdr1.PixelSpacing[0])*float(pre_hdr1.PixelSpacing[1])*float(abs(self.aff_mat[2,2])) #voxel size in mm^3

        voxsize_cm3 = voxsize_mm3/1000 #voxel size in cm^3
        print("voxel size in cm^3")
        print(voxsize_cm3)

        ftv_cm3 = ftv_voxels*voxsize_cm3 #ftv in cm^3
        ftv_cm3 = round(ftv_cm3,3)
        ftvcm3print = "Tumor volume: " + str(ftv_cm3) + " cc"
        self.ttlstr = ftvcm3print + "\n" + self.ttlstr_orig
        print("Tumor volume: " + str(ftv_voxels) + " voxels") #print ftv in cm^3 so that user can see it every time lesion segmentation is done
        print(ftvcm3print) #print ftv in cm^3 so that user can see it every time lesion segmentation is done

        self.heading = ftvcm3print


        #Manually update the lower left slice print info here so that this occurs as soon as segmentation is overlayed
        #and doesn't force user to slice scroll in order to show new FTV value
        #Update red and yellow slices with new text
        viewax = slicer.app.layoutManager().sliceWidget('Red').sliceView()
        viewax.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.LowerLeft,self.ttlstr) #Edit 10/30/2020: add exam details here
        viewax.forceRender()
            
        viewsag = slicer.app.layoutManager().sliceWidget('Yellow').sliceView()
        viewsag.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.LowerLeft,self.ttlstr) #Edit 10/30/2020: add exam details here
        viewsag.forceRender()

        

        #Update progress bar to say that FTV has been calculated
        progressBar.value = 80
        progressBar.labelText = "FTV calculated"
        slicer.app.processEvents()


        #add affine matrix IJKToRAS to numpy array so that it can be added to RGB image

    ##    volIJKToRASMat = vtk.vtkMatrix4x4()
    ##    inputVolume.GetIJKToRASMatrix(volIJKToRASMat)
    
    ##    aff_mat_RAS = np.zeros((4,4))
    ##    for i in range(4):
    ##      for j in range(4):
    ##        aff_mat_RAS[i,j] = volIJKToRASMat.GetElement(i,j)
        
        volRASToIJKMat = vtk.vtkMatrix4x4()
        inputVolume.GetRASToIJKMatrix(volRASToIJKMat)

        try:
          img_np_rgb = np.zeros((img_np.shape[0],img_np.shape[1],img_np.shape[2],3))
          img_np_rgb = img_np_rgb.astype('u1')

          #Use tumor mask and SER map to colorize input image
          img_np_rgb, nblue, pctblue, npurple, pctpurple, ngreen, pctgreen, nred, pctred, nyellow, pctyellow = ftv_map_gen_p2.serColorize3D(img_np_rgb,tumor_mask,voi_mask,self.ser,zs,zf,ys,yf,xs,xf)

          #Create strings containing the volume and % for each SER color
          bluestr = str(round(nblue*voxsize_cm3,3)) + " cc (" + str(round(pctblue,2)) + "%)"
          purplestr = str(round(npurple*voxsize_cm3,3)) + " cc (" + str(round(pctpurple,2)) + "%)"
          greenstr = str(round(ngreen*voxsize_cm3,3)) + " cc (" + str(round(pctgreen,2)) + "%)"
          redstr = str(round(nred*voxsize_cm3,3)) + " cc (" + str(round(pctred,2)) + "%)"
          yellowstr = str(round(nyellow*voxsize_cm3,3)) + " cc (" + str(round(pctyellow,2)) + "%)"
          totvolstr = str(ftv_cm3) + " cc"

          #Add these strings to module widget
          self.BlueValue.setText( bluestr )
          self.PurpleValue.setText( purplestr )
          self.GreenValue.setText( greenstr )
          self.RedValue.setText( redstr )
          self.YellowValue.setText( yellowstr )
          self.TotVolValue.setText( totvolstr )

          img_np_rgb = np.transpose(img_np_rgb,(2,1,0,3)) #must transpose to numpy z,y,x,(color) orientation
          #img_np_rgb[300:450,300:450,10:90,:] = [255,0,0] #test red region

          #Update progress bar to say SER colorized image has been created
          progressBar.value = 99
          progressBar.labelText = "SER colorized image created"
          slicer.app.processEvents()

          print("Tumor Percentage Distribution by SER Color")
          print("Blue:")
          print(pctblue)
          print("\nPurple:")
          print(pctpurple)
          print("\nGreen:")
          print(pctgreen)
          print("\nRed:")
          print(pctred)
          print("\nYellow:")
          print(pctyellow)


      ##    rgb_dtype = np.dtype([('R', 'u1'), ('G', 'u1'), ('B', 'u1')])
      ##    rgb_typed = img_np_rgb.view(rgb_dtype).reshape(shape_3d)
      
          #rgb_typed = np.transpose(rgb_typed,(1,0,2)) #apparently you need x and y dimension order to be switched to make orientation match with other images
          #rgb_typed = np.flip(rgb_typed,2) #apparently you also need to flip order of z-slices to make orientation match
      
      #try displaying rgb image directly instead of saving to .nii first
      #when I tried, I got the error: TypeError: Could not find a suitable VTK type for [('R', 'u1'), ('G', 'u1'), ('B', 'u1')]
      ##    segment_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")

          #Edit 7/10/2020: Andras Lasso said I should use vtkMRMLVectorVolumeNode to display SER colorized
          #lesion segmention directly without saving image first-- it worked!
          self.segment_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLVectorVolumeNode","SER Colorized Lesion Segmentation")
          slicer.util.updateVolumeFromArray(self.segment_node,img_np_rgb)
          self.segment_node.SetRASToIJKMatrix(volRASToIJKMat)
          slicer.util.setSliceViewerLayers(foreground=self.segment_node,foregroundOpacity = self.fopac)
          slicer.util.setSliceViewerLayers(background=inputVolume) #check if this helps you get rid of glitch of showing wrong background image when there are multiple visits
  ##        slicer.util.resetSliceViews()

          #Update progress bar to say lesion segmentation process is done
          progressBar.value = 100
          progressBar.labelText = "SER colorized image loaded to Slicer"
          slicer.app.processEvents()

          #Create copies of maps that can be reused for report creation
          self.sega = a
          self.segb = b
          self.segser = self.ser
          self.segtumormask = tumor_mask
          self.segvoimask = voi_mask

        except:
          slicer.util.confirmOkCancelDisplay("Only showing FTV value. Unable to display SER colorized map because input image is too large.","Error Displaying SER Colorization")


    
    ##    seg_savepath = qt.QFileDialog.getExistingDirectory(0,("Select save location for SER colorized image"))
    ##    rgb_img = nib.Nifti1Image(rgb_typed, aff_mat_RAS) #add IJKToRAS matrix to the rgb image to save
    ##    rgb_img_savename = seg_savepath + "\\" + "SERcolorized.nii"
    ##    nib.save(rgb_img,rgb_img_savename)
        #print("SER colorized image saved")

    
        #This does not work for 3D RGB image
    ##    segment_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    ##    slicer.util.updateVolumeFromArray(segment_node,img_np_rgb)
    ##    slicer.util.setSliceViewerLayers(background=segment_node)
    ##    slicer.util.resetSliceViews()

    #Edit 7/24/2020: If segmentLesion box is unchecked, remove vector volume
    else:
      slicer.util.setSliceViewerLayers(foreground=None,foregroundOpacity = self.fopac) #check if this helps you get rid of glitch after removing SER color when there are multiple visits
      #slicer.mrmlScene.RemoveNode(self.segment_node)



  def saveRegionToXML(self):
    #Edit 6/26/2020: Must make this avoid using LPS aff_mat and only use the image's RAS/IJK transform matrices


    #update datetime for filename every time you save a new xml file
    #start = datetime.datetime.now()
    #Format of runstarttime is like 20200101_12:30
    #self.runstarttime = str(start.year) + str(start.month) + str(start.day) + "_" + str(start.hour) + "_" + str(start.minute)  #Processing start time, that is included in output files or folders

    #Edit 7/6/2020: changing save file method so that path is automatically set to voi_lps_files and
    #user can specify name of xml file

    self.dflt_xml_svdir = self.exampath + "\\" + "voi_lps_files"

    #Make a new voi_lps_files folder for this exam if there isn't an existing one
    if(os.path.isdir(self.dflt_xml_svdir) == 0):
      os.mkdir(self.dflt_xml_svdir)

    dflt_xml_svnm = self.dflt_xml_svdir + "\\" + self.sitestr + "_" + self.idstr + "_" + self.nodevisstr + "_roiandomits.xml"

    savename = qt.QFileDialog.getSaveFileName(None, 'Save ROI XML File',
                                       dflt_xml_svnm,
                                       'XML (*.xml)')
    
    #savename = savepath + "\\" + "voi_omits_" + self.runstarttime + ".xml"

    inputVolume = self.inputSelector.currentNode()


    #read ROI and omit region info
    self.roi.GetXYZ(self.roicenter[0,:])
    
    self.roi.GetRadiusXYZ(self.roiradius[0,:])

    #convert ROI coordinates from RAS to IJK
##    self.roicenter[0,:], self.roiradius[0,:] = RASToIJKFunc(self.roicenter[0,:],self.roiradius[0,:],inputVolume)


    self.omit1.GetXYZ(self.omitcenters[0,:])
    self.omit2.GetXYZ(self.omitcenters[1,:])
    self.omit3.GetXYZ(self.omitcenters[2,:])
    self.omit4.GetXYZ(self.omitcenters[3,:])
    self.omit5.GetXYZ(self.omitcenters[4,:])

    self.omit1.GetRadiusXYZ(self.omitradii[0,:])
    self.omit2.GetRadiusXYZ(self.omitradii[1,:])
    self.omit3.GetRadiusXYZ(self.omitradii[2,:])
    self.omit4.GetRadiusXYZ(self.omitradii[3,:])
    self.omit5.GetRadiusXYZ(self.omitradii[4,:])

    #convert all omit coordinates from RAS to IJK
##    self.omitcenters[0,:], self.omitradii[0,:] = RASToIJKFunc(self.omitcenters[0,:],self.omitradii[0,:],inputVolume)
##    self.omitcenters[1,:], self.omitradii[1,:] = RASToIJKFunc(self.omitcenters[1,:],self.omitradii[1,:],inputVolume)
##    self.omitcenters[2,:], self.omitradii[2,:] = RASToIJKFunc(self.omitcenters[2,:],self.omitradii[2,:],inputVolume)
##    self.omitcenters[3,:], self.omitradii[3,:] = RASToIJKFunc(self.omitcenters[3,:],self.omitradii[3,:],inputVolume)
##    self.omitcenters[4,:], self.omitradii[4,:] = RASToIJKFunc(self.omitcenters[4,:],self.omitradii[4,:],inputVolume)


    #convert this info into voxel coordinates

    #get xs,xf,ys,yf,zs,zf,oxs,oxf,oys,oyf,ozs,ozf from center and radius variables
##    center = np.transpose(self.roicenter)
##    radius = np.transpose(self.roiradius)
##      
##    xs = int(center[0]) - int(radius[0])
##    xf = int(center[0]) + int(radius[0])
##
##    ys = int(center[1]) - int(radius[1])
##    yf = int(center[1]) + int(radius[1])
##       
##    zs = int(center[2]) - int(radius[2])
##    zf = int(center[2]) + int(radius[2])
##
##    #if no omits, all these variables are -1
##    if self.omitCount == 0:
##      oxs = np.array([-1])
##      oxf = np.array([-1])
##      oys = np.array([-1])
##      oyf = np.array([-1])
##      ozs = np.array([-1])
##      ozf = np.array([-1])
##    #otherwise, omit voxel variables should be made based on user defined omit regions only (using omitCount)
##    else:
##      oxs = self.omitcenters[0:self.omitCount,0] - self.omitradii[0:self.omitCount,0]
##      oxf = self.omitcenters[0:self.omitCount,0] + self.omitradii[0:self.omitCount,0]
##
##      oys = self.omitcenters[0:self.omitCount,1] - self.omitradii[0:self.omitCount,1]
##      oyf = self.omitcenters[:,1] + self.omitradii[:,1]
##       
##      ozs = self.omitcenters[0:self.omitCount,2] - self.omitradii[0:self.omitCount,2]
##      ozf = self.omitcenters[0:self.omitCount,2] + self.omitradii[0:self.omitCount,2]

    #create xml text and save to file
    if self.omitCount > 0:
      result = Write_to_xml.saveRASCoordsToLPSxml(self.roicenter,self.roiradius,self.omitcenters[0:self.omitCount,:],self.omitradii[0:self.omitCount,:],self.studydate)
    else:
      result = Write_to_xml.saveRASCoordsToLPSxml(self.roicenter,self.roiradius,np.array([-1]),np.array([-1]),self.studydate)

    file1 = open(savename,"w")
    file1.write(result)
    file1.close()
    
    print("ROI and omits saved to file")

    
#
# roiAndOmitsTestLogic
#

class DCE_TumorMapProcessLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  
  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True


  def run(self, gzipped, exampath, inputVolume,roiradius,roicenter,omitCount,omitradii,omitcenters,ijkToRASmat,tempres, fsort, nslice, manufacturer, studydate, dce_folders, aff_mat, aff_inv_mat, earlyPostContrastNum, latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss,sitestr,idstr,nodevisstr,window,level,pct,pethresh,minconnpix):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param outputVolume: thresholding result
    :param imageThreshold: values above/below this threshold will be set to 0
    :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
    :param showResult: show output volume in slice viewers
    """
    print("ROI and omit radii")
    print(roiradius)
    print(omitradii)

    print("ROI and omit centers")
    print(roicenter)
    print(omitcenters)

    #Edit 7/23/2020: No longer outputting ROI focused images
    

    #Function that applies voi mask to image and crops it to only include bounds of ROI
    def makeROIFocusedImg(img,voi_mask,xs,xf,ys,yf,zs,zf):
      img_roi = img*voi_mask
      img_roi = img_roi[xs:xf,ys:yf,zs:zf]
      return img_roi

    #Function that makes a node for the numpy image so that it can be viewed in a dropdown menu in Slicer window
    def makeImageNode(img,nodename):
      img = np.transpose(img,(2,1,0)) #numpy images that are displayed in Slicer must have dimension order z,y,x to have same orientation as NII and DICOM files
      img_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode",nodename) #add this image to dropdown node with name precontrast
      slicer.util.updateVolumeFromArray(img_node, img)
      return img_node
      

    logging.info('Processing started')
    progressBar = slicer.util.createProgressDialog()
    progressBar.value = 1
    progressBar.labelText = 'Processing started'
    slicer.app.processEvents()
    print("Processing started")


    #Edit 7/20/2020: If exam is gzipped, use gunzipped folder in exam folder for these functions
    if(gzipped == 1 and 'gunzipped' not in exampath):
      imagespath = exampath + "\\" + "gunzipped"
    else:
      imagespath = exampath

    #generate FTV processing outputs
    a,b,c,pe,ser,tumor_mask,voi_mask,zs,zf,ys,yf,xs,xf,pct,pre_thresh,pethresh,minconnpix = ftv_map_gen.makeFTVMaps(imagespath, manufacturer, dce_folders,roicenter,roiradius,omitcenters,omitradii, earlyPostContrastNum, latePostContrastNum,pct,pethresh,minconnpix)
      
    progressBar.value = 50
    progressBar.labelText = 'Computed PE, SER, and Tumor Mask'
    slicer.app.processEvents()
    print("Computed PE, SER, and Tumor Mask")


    #savename for figure
    #Edit 7/6/2020: changing save file method so that path is automatically set to voi_lps_files and
    #user can specify name of xml file

    #Define savepath for report and create new folder if necessary
    reportsavepath = exampath + "\\" + "SlicerReports"
    if(os.path.isdir(reportsavepath) == 0):
        os.mkdir(reportsavepath)
    
    dflt_report_svnm = reportsavepath + "\\" + sitestr + "_" + idstr + "_" + nodevisstr[:-1] + ".pdf"

    #prompt user to enter savename of report PDF
    savenamepdf = qt.QFileDialog.getSaveFileName(None, 'Save Slicer report',
                                       dflt_report_svnm,
                                       'PDF (*.pdf)')


    #generate report
    #Edit 6/30/2020: No more imgReversed
    ftv_plots.createPDFreport(gzipped,exampath,savenamepdf,tempres,fsort,manufacturer,dce_folders,nslice,earlyPostContrastNum,latePostContrastNum, earlydiffmm, earlydiffss, latediffmm, latediffss, a,b,ser,tumor_mask,voi_mask,xs,xf,ys,yf,zs,zf,omitCount,omitradii,omitcenters,pct,pre_thresh,pethresh,minconnpix,aff_mat,ijkToRASmat,nodevisstr,window,level,idstr)

    #After report is created, delete gunzipped folder
    gzip_gunzip_pyfuncs.deleteGunzipped(exampath)

    progressBar.value = 99
    progressBar.labelText = 'Report saved'
    slicer.app.processEvents()
    print("Report saved")

    #for purpose of testing out pyradiomics, save these 2 nrrd's
##    pe_roi = makeROIFocusedImg(pe,voi_mask,xs,xf,ys,yf,zs,zf)
##    pe_roi_savename = exampath + "\\" + "PE_roi.nrrd"
##    nrrd.write(pe_roi_savename,pe_roi)
##
##    tumor_mask_roi = makeROIFocusedImg(tumor_mask,voi_mask,xs,xf,ys,yf,zs,zf)
##    tumor_mask_roi_savename = exampath + "\\" + "tumormask_roi.nrrd"
##    nrrd.write(tumor_mask_roi_savename,tumor_mask_roi)



    #make cropped pre-contrast and save it to node
##    a_roi = makeROIFocusedImg(a,voi_mask,xs,xf,ys,yf,zs,zf)
##    a_node = makeImageNode(a_roi,"Precontrast ROI")
##
##    progressBar.value = 50
##    progressBar.labelText = 'Pre-contrast ROI loaded to Slicer'
##    slicer.app.processEvents()
##    print("Pre-contrast ROI loaded to Slicer")
##
##
##    #make cropped PE and save it to node
##    pe_roi = makeROIFocusedImg(pe,voi_mask,xs,xf,ys,yf,zs,zf)
##    pe_node = makeImageNode(pe_roi,"PE ROI")
##
##    progressBar.value = 75
##    progressBar.labelText = 'PE ROI loaded to Slicer'
##    slicer.app.processEvents()
##    print("PE ROI loaded to Slicer")
##
##
##    #make cropped tumor mask and save it to node
##    tumor_mask_roi = makeROIFocusedImg(tumor_mask,voi_mask,xs,xf,ys,yf,zs,zf)
##    tumor_mask_node = makeImageNode(tumor_mask_roi,"Tumor Mask ROI")
##
##    progressBar.value = 99
##    progressBar.labelText = 'Tumor Mask ROI loaded to Slicer'
##    slicer.app.processEvents()
##    print("Tumor Mask ROI loaded to Slicer")

    time.sleep(1) #pause for 1 second at progress bar 99%

    progressBar.value = 100
    progressBar.labelText = 'Processing completed'
    slicer.app.processEvents()


    logging.info('Processing completed')




#
# roiAndOmitsTestTest
#

class DCE_TumorMapProcessTest(ScriptedLoadableModuleTest):
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
    self.test_roiAndOmitsTest1()

  def test_DCE_TumorMapProcess1(self):
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
    logic = roiAndOmitsTestLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
