find_program(qt_xmlpatterns_executable
  NAMES xmlpatterns-qt5 xmlpatterns
  DOC "Path to xmlpatterns")
mark_as_advanced(qt_xmlpatterns_executable)

paraview_add_plugin(pv_zoltan
  VERSION "1.0"
  MODULES pvZoltanFilters
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../Filters/vtk.module"
  XML_DOCUMENTATION OFF)
