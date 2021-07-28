cmake_minimum_required(VERSION 3.13.4)

project(Breast_DCEMRI_FTV)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/rnadkarni2/SlicerBreast_DCEMRI_FTV")
set(EXTENSION_CATEGORY "FTV Segmentation")
set(EXTENSION_CONTRIBUTORS "Rohan Nadkarni (UCSF Breast Imaging Research Group), David Newitt (UCSF Breast Imaging Research Group)")
set(EXTENSION_DESCRIPTION "This extension contains both of the modules needed to compute functional tumor volume (FTV) from a dynamic contrast-enhanced magnetic resonance imaging (DCE-MRI) exam of a breast cancer patient. The current version of this extension is intended for use on bilateral breast MR exams with axial slices.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/rnadkarni2/SlicerBreast_DCEMRI_FTV/master/Breast_DCEMRI_FTV.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/rnadkarni2/SlicerBreast_DCEMRI_FTV/master/Module1Screenshot.png https://raw.githubusercontent.com/rnadkarni2/SlicerBreast_DCEMRI_FTV/master/Module2Screenshot.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(DCE_IDandPhaseSelect)
add_subdirectory(DCE_TumorMapProcess)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
