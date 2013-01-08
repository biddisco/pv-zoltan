project("pv_zoltan")

#--------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
endif(COMMAND cmake_policy)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF (ParaView_SOURCE_DIR)
  IF(PARAVIEW_BUILD_QT_GUI)
    INCLUDE(${QT_USE_FILE})
  ENDIF(PARAVIEW_BUILD_QT_GUI)
  INCLUDE_DIRECTORIES(
    ${HDF5_INCLUDE_DIR}
    ${PARAVIEW_INCLUDE_DIRS}
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

#--------------------------------------------------
# Set project include directories 
# pv-common is assumed to be checkout out too
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
)

# --------------------------------------------------
# Trilinos support 
# --------------------------------------------------
SET(TRILINOS_ENABLED 0)
IF (NOT PARAVIEW_USE_MPI)
  SET(USE_TRILINOS OFF CACHE BOOL "Use Trilinos library" FORCE)
  message(FATAL_ERROR "pv-zoltan cannot be built against a non MPI enabled ParaView")
ENDIF (NOT PARAVIEW_USE_MPI)

OPTION(USE_SYSTEM_TRILINOS "Use System Trilinos library" OFF)
IF(USE_SYSTEM_TRILINOS)
  FIND_PACKAGE(Trilinos)
ELSE(USE_SYSTEM_TRILINOS)
  SET(TPL_ENABLE_MPI ON CACHE BOOL "Do not change")
  SET(Trilinos_ENABLE_TESTS OFF CACHE BOOL "Do not change")
  SET(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "Do not change")
  SET(Trilinos_ENABLE_Zoltan ON CACHE BOOL "Do not change")
  SET(Trilinos_ENABLE_Fortran OFF CACHE BOOL "Do not change")
  IF (WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include")
  ENDIF(WIN32)
  ADD_DEFINITIONS(-DREALLOC_BUG)
  ADD_SUBDIRECTORY(trilinos)
  SET(Trilinos_FOUND 1)
  SET(Trilinos_DEPS_HTML_OUTPUT_FILE "" CACHE STRING "Leave empty" FORCE)
ENDIF(USE_SYSTEM_TRILINOS)

IF(Trilinos_FOUND)
  INCLUDE_DIRECTORIES("${Trilinos_INCLUDE_DIRS}")
  LINK_DIRECTORIES("${Trilinos_LIBRARY_DIRS}")
  # MESSAGE("Trilinos Include Dirs ${Trilinos_INCLUDE_DIRS} and Trilinos Link Dirs ${Trilinos_LIBRARY_DIRS} ")
  #
  # Add source which depends on trilinos to our plugin
  #
  SET(TRILINOS_CXX 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanV1PartitionFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMeshPartitionFilter.cxx
  )
  SET(TRILINOS_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanV1PartitionFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMeshPartitionFilter.xml
  )
  SET(TRILINOS_LIBS zoltan)
  SET(TRILINOS_ENABLED 1)
  ADD_DEFINITIONS(-DPLUGIN_USE_TRILINOS)
ELSE(Trilinos_FOUND)
  message(FATAL_ERROR "pv-zoltan cannot be built unless Trilinos is correctly found")
ENDIF(Trilinos_FOUND)

# --------------------------------------------------
# Plugin config
# --------------------------------------------------
set(XML_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")
set(SRC_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")

#--------------------------------------------------
# Define plugin sources
#--------------------------------------------------
SET(PLUGIN_NAME pv_zoltan)
ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 

  SERVER_MANAGER_XML     
    ${TRILINOS_XML}
    ${XML_PATH}/vtkPartitionOutline.xml

  SERVER_MANAGER_SOURCES
    ${TRILINOS_CXX}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkBoundsExtentTranslator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPartitionOutline.cxx

  SERVER_SOURCES

  GUI_INTERFACES 
    ${ZOLTAN_IFACE}

  GUI_SOURCES
    ${ZOLTAN_IFACE_SRCS}
    ${ZOLTAN_MOC_SRCS}
    ${ZOLTAN_UI_SOURCES}

  GUI_RESOURCE_FILES
#    ${XML_PATH}/pv_ZOLTAN_Readers.xml
#    ${XML_PATH}/pv_ZOLTAN_Sources.xml
    ${XML_PATH}/pv_zoltan_Filters.xml
)

#--------------------------------------------------
# setup plugin linking
#--------------------------------------------------
TARGET_LINK_LIBRARIES(${PLUGIN_NAME} 
#        vtkRenderingFreeTypeFontConfig
    ${ZOLTAN_SDK_LIBRARIES}
    ${TRILINOS_LIBS}
    ${HDF5_C_LIBRARY} 
    ${HDF5_CXX_LIBRARY} 
    ${Boost_SERIALIZATION_LIBRARY} 
    ${Boost_THREAD_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
)

#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(testing)
ENDIF (BUILD_TESTING)  


