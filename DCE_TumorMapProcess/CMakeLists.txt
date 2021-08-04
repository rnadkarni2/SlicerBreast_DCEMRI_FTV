#-----------------------------------------------------------------------------
set(MODULE_NAME DCE_TumorMapProcess)

add_subdirectory(Breast_DCEMRI_FTV_plugins2)
#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  $Breast_DCEMRI_FTV_plugins2/read_DCE_images_to_numpy.py
  $Breast_DCEMRI_FTV_plugins2/multivolume_folder_sort.py
  $Breast_DCEMRI_FTV_plugins2/ident_gzipped_exam.py
  $Breast_DCEMRI_FTV_plugins2/gzip_gunzip_pyfuncs.py
  $Breast_DCEMRI_FTV_plugins2/ftv_plots.py
  $Breast_DCEMRI_FTV_plugins2/ftv_map_gen_p2.py
  $Breast_DCEMRI_FTV_plugins2/ftv_map_gen.py
  $Breast_DCEMRI_FTV_plugins2/create2DimgAllFunctions.py
  $Breast_DCEMRI_FTV_plugins2/compute_lps_to_rcs.py
  $Breast_DCEMRI_FTV_plugins2/Write_to_xml.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/UI/${MODULE_NAME}.ui
  $ser_colorbar.PNG
  $Breast_DCEMRI_FTV_plugins2/3DSlicerLogo.PNG
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
