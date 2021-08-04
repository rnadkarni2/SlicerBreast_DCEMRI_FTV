#-----------------------------------------------------------------------------
set(MODULE_NAME DCE_IDandPhaseSelect)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  Breast_DCEMRI_FTV_plugins1/__init__.py
  Breast_DCEMRI_FTV_plugins1/timing_info_class_all_manufacturer.py
  Breast_DCEMRI_FTV_plugins1/read_DCE_images_to_numpy.py
  Breast_DCEMRI_FTV_plugins1/multivolume_folder_sort.py
  Breast_DCEMRI_FTV_plugins1/ident_gzipped_exam.py
  Breast_DCEMRI_FTV_plugins1/gzip_gunzip_pyfuncs.py
  Breast_DCEMRI_FTV_plugins1/compute_lps_to_rcs.py
  Breast_DCEMRI_FTV_plugins1/chooseEarly_LateByManufacturer.py
  Breast_DCEMRI_FTV_plugins1/Get_header_info_all_manufacturer.py
  Breast_DCEMRI_FTV_plugins1/Exam_Ident_and_timing.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/UI/${MODULE_NAME}.ui
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
