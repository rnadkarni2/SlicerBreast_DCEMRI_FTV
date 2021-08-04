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
#File that is called to write new user-selected VOI and omit regions to xml file


def saveNewRegionsToXML(path,aff_mat,studydate, xsnew, xfnew, ysnew, yfnew, zsnew, zfnew, oxsnew, oxfnew, oysnew, oyfnew, ozsnew, ozfnew):
    import yattag
    from yattag import Doc, indent
    import numpy as np


    #Source: https://stackoverflow.com/questions/3605680/creating-a-simple-xml-file-using-python

    #Allow user to adjust VOI for selected path
    #and output adjusted VOI parameters in xml file format

    def voxelsToXmlInput(aff_mat, xsnew,xfnew,ysnew,yfnew,zsnew,zfnew):
        #----------------------------------------------------------------------------------------------------------------------------------------------------------
        #3. Convert new VOI from RCS (row, column, slice) coordinates to LPS coordinates
        start_rcshomog = np.zeros((4,1))
        start_rcshomog[0,0] = ysnew
        start_rcshomog[1,0] = xsnew
        start_rcshomog[2,0] = zsnew
        start_rcshomog[3,0] = 1
        start_lpshomog = np.matmul(aff_mat,start_rcshomog)

        final_rcshomog = np.zeros((4,1))
        final_rcshomog[0,0] = yfnew
        final_rcshomog[1,0] = xfnew
        final_rcshomog[2,0] = zfnew
        final_rcshomog[3,0] = 1
        final_lpshomog = np.matmul(aff_mat,final_rcshomog)

        centerlps = (start_lpshomog+final_lpshomog)/2
        centerlps = centerlps[0:3]
        centerlps = np.transpose(centerlps)

        widthhalflencoords = np.zeros((1,3))
        widthhalflencoords.astype('float64')
        widthhalflencoords[0,0] = start_lpshomog[0,0]-float(centerlps[0,0])

        heighthalflencoords = np.zeros((1,3))
        heighthalflencoords.astype('float64')
        heighthalflencoords[0,1] = start_lpshomog[1,0]-float(centerlps[0,1])

        depthhalflencoords = np.zeros((1,3))
        depthhalflencoords.astype('float64')
        depthhalflencoords[0,2] = start_lpshomog[2,0]-float(centerlps[0,2])

        return centerlps, widthhalflencoords, heighthalflencoords, depthhalflencoords

    #Convert to xml input format for VOI
    centerlps, widthhalflencoords, heighthalflencoords, depthhalflencoords = voxelsToXmlInput(aff_mat,xsnew,xfnew,ysnew,yfnew,zsnew,zfnew)



    #Define xml file format
    doc, tag, text = Doc().tagtext()
    omit = 1
    doc.asis('<?xml version="1.0" encoding="utf-8"?>')
    text('\n')
    with tag('AegisXmlObject'):
        doc.stag('VersionInfo ')
        with tag('RoiSavedState'):
            with tag('StudyDate'):
                text(studydate)
            with tag('RoiInfo'):
                with tag('Item-CubicRoiInfo3D'):
                    doc.stag('Center',X=str(centerlps[0,0]),Y=str(centerlps[0,1]),Z=str(centerlps[0,2]))
                    doc.stag('WidthHalfLength',X=str(widthhalflencoords[0,0]),Y=str(widthhalflencoords[0,1]),Z=str(widthhalflencoords[0,2]))
                    doc.stag('HeightHalfLength',X=str(heighthalflencoords[0,0]),Y=str(heighthalflencoords[0,1]),Z=str(heighthalflencoords[0,2]))
                    doc.stag('DepthHalfLength',X=str(depthhalflencoords[0,0]),Y=str(depthhalflencoords[0,1]),Z=str(depthhalflencoords[0,2]))
                    with tag('OmitRegion'):
                        text('False')
                #adding an 'Item-CubicRoiInfo3D' for any omit regions that may exist
                if(oxsnew[0]!=-1):
                    ocenterlps = np.zeros((len(oxsnew),3))
                    owidthhalflencoords = np.zeros((len(oxsnew),3))
                    oheighthalflencoords = np.zeros((len(oxsnew),3))
                    odepthhalflencoords = np.zeros((len(oxsnew),3))

                    for i in range(len(oxsnew)):
                        ocenterlps[i,:], owidthhalflencoords[i,:], oheighthalflencoords[i,:], odepthhalflencoords[i,:] = voxelsToXmlInput(aff_mat,oxsnew[i],oxfnew[i],oysnew[i],oyfnew[i],ozsnew[i],ozfnew[i])
                        with tag('Item-CubicRoiInfo3D'):
                            doc.stag('Center',X=str(ocenterlps[i,0]),Y=str(ocenterlps[i,1]),Z=str(ocenterlps[i,2]))
                            doc.stag('WidthHalfLength',X=str(owidthhalflencoords[i,0]),Y=str(owidthhalflencoords[i,1]),Z=str(owidthhalflencoords[i,2]))
                            doc.stag('HeightHalfLength',X=str(oheighthalflencoords[i,0]),Y=str(oheighthalflencoords[i,1]),Z=str(oheighthalflencoords[i,2]))
                            doc.stag('DepthHalfLength',X=str(odepthhalflencoords[i,0]),Y=str(odepthhalflencoords[i,1]),Z=str(odepthhalflencoords[i,2]))
                            with tag('OmitRegion'):
                                text('True')


    result = indent(
        doc.getvalue(),
        indentation = ' '*2,
        newline = '\n'
    )


    return result


def saveRASCoordsToLPSxml(roicenter,roiradius,ocenters,oradii,studydate):
    import yattag
    from yattag import Doc, indent
    import numpy as np


    try:
        if ocenters[0,0] != -1:
            omits = 1
        else:
            omits = 0
    except:
        if ocenters[0] == -1:
            omits = 0
        else:
            omits = 1

    #Define xml file format
    doc, tag, text = Doc().tagtext()
    omit = 1
    doc.asis('<?xml version="1.0" encoding="utf-8"?>')
    text('\n')
    with tag('AegisXmlObject'):
        doc.stag('VersionInfo ')
        with tag('RoiSavedState'):
            with tag('StudyDate'):
                text(studydate)
            with tag('RoiInfo'):
                with tag('Item-CubicRoiInfo3D'):
                    doc.stag('Center',X=str(-roicenter[0,0]),Y=str(-roicenter[0,1]),Z=str(roicenter[0,2]))
                    doc.stag('WidthHalfLength',X=str(roiradius[0,0]),Y=str(0),Z=str(0))
                    doc.stag('HeightHalfLength',X=str(0),Y=str(roiradius[0,1]),Z=str(0))
                    doc.stag('DepthHalfLength',X=str(0),Y=str(0),Z=str(roiradius[0,2]))
                    with tag('OmitRegion'):
                        text('False')
                #adding an 'Item-CubicRoiInfo3D' for any omit regions that may exist
                if(omits == 1):

                    for i in range(ocenters.shape[0]):
                        with tag('Item-CubicRoiInfo3D'):
                            doc.stag('Center',X=str(-ocenters[i,0]),Y=str(-ocenters[i,1]),Z=str(ocenters[i,2]))
                            doc.stag('WidthHalfLength',X=str(oradii[i,0]),Y=str(0),Z=str(0))
                            doc.stag('HeightHalfLength',X=str(0),Y=str(oradii[i,1]),Z=str(0))
                            doc.stag('DepthHalfLength',X=str(0),Y=str(0),Z=str(oradii[i,2]))
                            with tag('OmitRegion'):
                                text('True')


    result = indent(
        doc.getvalue(),
        indentation = ' '*2,
        newline = '\n'
    )


    return result

