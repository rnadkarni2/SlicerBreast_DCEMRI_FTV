<h1>Breast DCE-MRI FTV extension</h1>

<p>Users of this extension must abide by the Slicer license. A copy of this license is included
in this repository in the file <strong>license.txt</strong>. The <strong>copyright.txt</strong> 
file in this repository includes the copyright statement written on all Python code files in this extension.</p>

<h2>Purpose of Extension</h2>

<p>This 3D Slicer extension allows the user to compute functional tumor volume (FTV) from a 
dynamic contrast-enhanced magnetic resonance imaging (DCE-MRI) exam of a breast 
cancer patient. FTV is a quantitative measure of tumor burden. It is calculated through a 
segmentation process developed by the UCSF Breast Imaging Research Group that accounts 
for pre-contrast enhancement and post-contrast kinetics.</p>

<h2>Background</h2>

<p>FTV from DCE-MRI has already been established as a reliable biomarker in 
breast cancer patients undergoing neoadjuvant chemotherapy (NAC). Research from the
ACRIN 6657 / I-SPY 1 TRIAL showed that FTV is a stronger predictor of  treatment response 
than clinical assessment. Due to this promising result, the ongoing I-SPY 2 TRIAL is using FTV 
to guide adaptive randomization of therapy for patients.</p>

<p>Despite the predictive power of the FTV segmentation method, access to it has been limited.
The most recent implementation of this method prior to this 3D Slicer extension, a module on the
Sentinelle AEGIS platform, can only be used by sites participating in the I-SPY 2 TRIAL due to
restrictions of the FDA IDE approval. Therefore, the Breast DCE-MRI FTV extension
has been developed and published so that the larger breast MR imaging community can use the
FTV segmentation method by simply downloading 3D Slicer and installing this extension
through the Extension Manager.</p>


<h2>Summary of Modules</h2>

<p><strong>This extension is only compatible with MR exams that use bilateral
images with axial slices. It has been used to process numerous MR exams on
a Windows computer and has been shown to work on a Mac computer as well.</strong> </p>

<p>This extension contains 2 modules. After you install the extension, both of these modules
can be found under the <strong>FTV Segmentation</strong> category in the modules list.</p>

<p>The first module is called <strong>Module 1: Load DCE Images</strong>. In this module, the user selects the MR exam to process,
preferences for early and late phase timing, and method of DCE series identification
(automatic or manual). The module will then load the pre-contrast, early post-contrast,
and late post-contrast phases to the Slicer window. </p>

![Screenshot from Module 1](https://github.com/rnadkarni2/SlicerBreast_DCEMRI_FTV/blob/master/Module1Screenshot.png)

<p>The second module is called <strong>Module 2: Compute FTV in ROI</strong>. In this module, the user first selects the tumor region of interest
(ROI) and also has the option to select regions within the ROI to exclude from the segmentation (omit regions).
The phases loaded by Module 1 as well as subtraction images and MIPs derived from
them may be used to guide selection of ROI and omit regions. Furthermore, the user has the option
to either create ROI and omit regions by dragging bounding boxes to the appropriate position and size 
or import them from an existing .xml file. If the user creates a new
ROI and omit regions, they can be saved to an xml file for future use.</p>

<p>After ROI and omit region selections, the user can begin the segmentation process.
If necessary, the user may change the segmentation thresholds from default values.
The user can segment the lesion and overlay a tumor segmentation that is color-coded
by signal enhancement ratio (SER) values onto the axial and sagittal views.
In addition, the user may generate a report that includes several axial and sagittal images from the exam
as well as relevant information such as FTV value, ROI boundaries, segmentation thresholds used,
and details about the exam selected.</p>

![Screenshot from Module 2](https://github.com/rnadkarni2/SlicerBreast_DCEMRI_FTV/blob/master/Module2Screenshot.png)

<h2>Publication</h2>

<p>A publication about this 3D Slicer extension called <strong>Validation of the Open-Source 3D Slicer Extension Breast DCE-MRI Functional Tumor Volume as a Tool to Assess Treatment Response</strong> 
will be submitted to the journal <strong>Tomography</strong> before the end of 2021.</p>

<h2>Additional Resources</h2>

For complete, step-by-step user instructions on how to use this extension, please read the document <br>
**UserInstructionsFor_Breast_DCEMRI_FTV.pdf**.

<p>If you have any questions, you can contact David Newitt of the UCSF Breast Imaging Research Group (david.newitt@ucsf.edu)</p>

<h2>References</h2>
<p>1. Hylton NM, Blume JD, Bernreuter WK, Pisano ED, Rosen MA, Morris EA, Weatherall PT, Lehman CD, Newstead GM, Polin S. Locally advanced breast cancer: MR imaging for prediction of response to neoadjuvant chemotherapy—results from ACRIN 6657/I-SPY TRIAL. Radiology. 2012;263:663–672.</p>
<p>2. Hylton NM. Vascularity assessment of breast lesions with gadolinium-enhanced MR imaging. Magn Reson Imaging Clin N Am. 1999;7:411–420. x. </p>