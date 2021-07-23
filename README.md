<h1>Breast DCE-MRI FTV extension</h1>

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
restrictions of the FDA IDE approval. Therefore, the Breast DCE-MRI Tumor Segment extension
has been developed and published so that the larger breast MR imaging community can use the
FTV segmentation method by simply downloading 3D Slicer and installing this extension
through the Extension Manager.</p>


<h2>Summary of Modules</h2>

<p>This extension contains 2 modules. After you install the extension, both of these modules
can be found under the <strong>FTV Segmentation</strong> category in the modules list.</p>

<p>The first module is called <strong>Module 1: Load DCE Images</strong>. In this module, the user selects the MR exam to process,
preferences for early and late phase timing, and method of DCE series identification
(automatic or manual). The module will then load the pre-contrast, early post-contrast,
and late post-contrast phases to the Slicer window. <strong>This extension is only compatible with MR exams that use bilateral
images with axial slices.</strong> </p>

![Screenshot from Module 1](https://github.com/rnadkarni2/SlicerBreast_DCEMRI_TumorSegment/blob/master/Module1Screenshot.png)

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

![Screenshot from Module 2](https://github.com/rnadkarni2/SlicerBreast_DCEMRI_TumorSegment/blob/master/Module2Screenshot.png)


<h2>Additional Resources</h2>

For complete, step-by-step user instructions on how to use this extension, please read the document <br>
**UserInstructionsFor_Breast_DCEMRI_FTV.pdf**.

For more information about the FTV segmentation method, please read the publication <br>
***Vascularity assessment of breast lesions with gadolinium-enhanced MR imaging***.