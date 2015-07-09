project("pv_zoltan")

#--------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#--------------------------------------------------
# cmake policies
# see http://www.cmake.org/cmake/help/v2.8.12/cmake.html#section_Policies
#--------------------------------------------------
foreach(p
    CMP0003 
    CMP0017 
    CMP0023 
    CMP0042
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

add_definitions(-DCMAKE_IS_BROKEN)

#-----------------------------------------------------------------------------
# Testing support
#-----------------------------------------------------------------------------
include(CTest)

#-----------------------------------------------------------------------------
# Basic naming/config stuff here
#-----------------------------------------------------------------------------
set(pv_zoltan_VERSION_MAJOR   0)
set(pv_zoltan_VERSION_MINOR   9)
set(pv_zoltan_VERSION_RELEASE 0)

set (pv_zoltan                "pv_zoltan")
set (pv_zoltan_NAME           "pv_zoltan")
set (pv_zoltan_VERSION        "${pv_zoltan_VERSION_MAJOR}.${pv_zoltan_VERSION_MINOR}.${pv_zoltan_VERSION_RELEASE}")
set (pv_zoltan_VERSION_STRING "${pv_zoltan_VERSION}")
set (pv_zoltan_SOVERSION      "${pv_zoltan_SOVERS_MAJOR}.${pv_zoltan_SOVERS_RELEASE}.${pv_zoltan_SOVERS_MINOR}")
set (pv_zoltan_STRING         "${pv_zoltan_NAME} ${pv_zoltan_VERSION_STRING}")
set (pv_zoltan_TARNAME        "${pv_zoltan}")
set (pv_zoltan_URL            "http://www.cscs.ch")
set (pv_zoltan_BUGREPORT      "help@cscs.ch")

#--------------------------------------------------
# Location of cmake files
#--------------------------------------------------
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
)

#------------------------------------------------------
# MPI support
#------------------------------------------------------
IF (PARAVIEW_USE_MPI)
  FIND_PACKAGE(MPI)   
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ELSE (PARAVIEW_USE_MPI)
  SET(USE_TRILINOS OFF CACHE BOOL "Use Trilinos library" FORCE)
  message(FATAL_ERROR "pv-zoltan requires an MPI enabled ParaView")
ENDIF (PARAVIEW_USE_MPI)

# --------------------------------------------------
# We must build sub libraries with -fPIC
# --------------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#--------------------------------------------------
# GitExternal config
#--------------------------------------------------
set(GIT_EXTERNAL_DISABLE_UPDATE ON)
set(GIT_EXTERNAL_USER_FORK "" CACHE STRING "Disabled")

#--------------------------------------------------
# load script for checking out projects from git
#--------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/GitExternal.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SubProject.cmake)

#--------------------------------------------------
# checkout reduced version of trilinos 
#--------------------------------------------------
Git_External(
  ${CMAKE_CURRENT_SOURCE_DIR}/trilinos
  https://github.com/biddisco/trilinos
  zoltan_only 
  NO_UPDATE
  VERBOSE
)

# --------------------------------------------------
# Trilinos support 
# --------------------------------------------------
SET(TRILINOS_ENABLED 0)
OPTION(USE_SYSTEM_TRILINOS "Use System Trilinos library" OFF)
IF(USE_SYSTEM_TRILINOS)
  FIND_PACKAGE(Trilinos)
ELSE(USE_SYSTEM_TRILINOS)
  set(Trilinos_ASSERT_MISSING_PACKAGES      OFF CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Eptra                 ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Tpetra                 ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Xpetra                 ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Teuchos                ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Kokkos                 ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Zoltan                 ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Zoltan2                ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Fortran               OFF CACHE BOOL "Do not change")
  set(TPL_ENABLE_MPI                         ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_TESTS                 OFF CACHE BOOL "Do not change")
  IF (WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include")
  ENDIF(WIN32)
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
    ${Module_SRCS}
    ${Module_CUDA_C_SRCS}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPKdTree2.cxx

  GUI_INTERFACES 
    ${ZOLTAN_IFACE}

  GUI_SOURCES
    ${ZOLTAN_IFACE_SRCS}
    ${ZOLTAN_MOC_SRCS}
    ${ZOLTAN_UI_SOURCES}

  GUI_RESOURCE_FILES
    ${XML_PATH}/pv_zoltan_Filters.xml
)

#--------------------------------------------------
# setup plugin linking
#--------------------------------------------------
TARGET_LINK_LIBRARIES(${PLUGIN_NAME} 
  PUBLIC
    vtksys
    ${ZOLTAN_SDK_LIBRARIES}
    ${TRILINOS_LIBS}
)

#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(testing)
ENDIF (BUILD_TESTING)  

#--------------------------------------------------------
# Create the configuration header.
#--------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/module_configure.h.in
               ${PROJECT_BINARY_DIR}/pv_zoltan_configure.h
               @ONLY IMMEDIATE)

#-----------------------------------------------------------------------------
# Configure the project-config.cmake file for the build directory
#-----------------------------------------------------------------------------
set (${PROJECT_NAME}_INCLUDES_BUILD_TIME
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_BINARY_DIR}/trilinos/packages/zoltan/src"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan/src/include"
)
if(WIN32)
  set(${PROJECT_NAME}_INCLUDES_BUILD_TIME 
    ${PROJECT_SOURCE_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include
    ${${PROJECT_NAME}_INCLUDES_BUILD_TIME}      
  )
endif()


return()

include(CMakeInstallation)

#--------------------------------------------------
# Install
#--------------------------------------------------
SET(INSTALL_PATH 
  "${CMAKE_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}"
)

INSTALL(
  TARGETS 
    ${PLUGIN_NAME}
  DESTINATION ${INSTALL_PATH}
)
