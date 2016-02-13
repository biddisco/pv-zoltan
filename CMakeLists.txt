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

#--------------------------------------------------
# Setup Output directories
#--------------------------------------------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
        ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Executables."
        )
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
        ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Libraries"
        )
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all static libraries."
        )

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
  zoltan_subset
  NO_UPDATE
  VERBOSE
)

#--------------------------------------------------
# option to use old zoltan 1 version
#--------------------------------------------------
option(PV_ZOLTAN_USE_ZOLTAN_1 "Use version 1 of zoltan (legacy)" OFF)
if(PV_ZOLTAN_USE_ZOLTAN_1)
  add_definitions(-DVTK_ZOLTAN1_PARTITION_FILTER)
else()
  add_definitions(-DVTK_ZOLTAN2_PARTITION_FILTER)
endif()

#--------------------------------------------------
# CLAPACK
#--------------------------------------------------
if (WIN32 AND NOT PV_ZOLTAN_USE_ZOLTAN_1)
  include(${CMAKE_CURRENT_SOURCE_DIR}/trilinos/cmake/tribits/core/utils/AdvancedSet.cmake)
  add_subdirectory(clapack-3.2.1-CMAKE)
  set(CLAPACK_DIR ${PROJECT_BINARY_DIR}/clapack-3.2.1-CMAKE CACHE STRING "Do not change" FORCE)
  FIND_PACKAGE(CLAPACK 3.2.1 NO_MODULE)
    IF (CLAPACK_FOUND)
    ADVANCED_SET(TPL_BLAS_LIBRARIES blas
      CACHE FILEPATH "Set from MSVC CLAPACK specialization")
  ENDIF()
  set(BLAS_LIBRARY_DIRS "${PROJECT_BINARY_DIR}/bin")
  set(LAPACK_LIBRARY_DIRS "${PROJECT_BINARY_DIR}/bin")
  set(TPL_BLAS_LIBRARIES "blas")
  set(TPL_LAPACK_LIBRARIES "lapack")
#  install(EXPORT "blas" )
endif()

# --------------------------------------------------
# Trilinos support 
# --------------------------------------------------
OPTION(USE_SYSTEM_TRILINOS "Use System Trilinos library" OFF)
IF(USE_SYSTEM_TRILINOS)
  FIND_PACKAGE(Trilinos)
ELSE(USE_SYSTEM_TRILINOS)
  set(Trilinos_ASSERT_MISSING_PACKAGES      OFF CACHE BOOL "Do not change")

  set(Trilinos_ENABLE_Zoltan                 ON CACHE BOOL "Do not change" FORCE)
  if (PV_ZOLTAN_USE_ZOLTAN_1)
    set(Trilinos_ENABLE_Tpetra                 OFF CACHE BOOL "Do not change" FORCE)
    set(Trilinos_ENABLE_Zoltan2                OFF CACHE BOOL "Do not change" FORCE) 
    set(TPL_ENABLE_BLAS                        OFF CACHE BOOL "Do not change" FORCE)
    set(TPL_ENABLE_LAPACK                      OFF CACHE BOOL "Do not change" FORCE)
  else()
    set(Trilinos_ENABLE_Tpetra                 ON CACHE BOOL "Do not change" FORCE)
    set(Trilinos_ENABLE_Zoltan2                ON CACHE BOOL "Do not change" FORCE) 
    set(TPL_ENABLE_BLAS                        ON CACHE BOOL "Do not change" FORCE)
    set(TPL_ENABLE_LAPACK                      ON CACHE BOOL "Do not change" FORCE)
  endif()

  set(TPL_ENABLE_MPI                         ON CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_Fortran               OFF CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "Do not change")
  set(Trilinos_ENABLE_TESTS                 OFF CACHE BOOL "Do not change")
  set(Trilinos_REQUIRES_PYTHON              OFF CACHE BOOL "Do not change")
  set(Trilinos_USES_PYTHON                  OFF CACHE BOOL "Do not change")

  IF (WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include")
  ENDIF(WIN32)
  set(BUILD_SHARED_LIBS OFF)
  ADD_SUBDIRECTORY(trilinos)
  SET(Trilinos_FOUND 1)
  SET(Trilinos_DEPS_HTML_OUTPUT_FILE "" CACHE STRING "Leave empty" FORCE)
ENDIF(USE_SYSTEM_TRILINOS)


IF(NOT Trilinos_FOUND)
  message(FATAL_ERROR "pv-zoltan cannot be built unless Trilinos is correctly found")
ENDIF(Trilinos_FOUND)

# --------------------------------------------------
# Include dirs for plugin
# --------------------------------------------------
INCLUDE_DIRECTORIES("${Trilinos_INCLUDE_DIRS}")
LINK_DIRECTORIES("${Trilinos_LIBRARY_DIRS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# MESSAGE("Trilinos Include Dirs ${Trilinos_INCLUDE_DIRS} and Trilinos Link Dirs ${Trilinos_LIBRARY_DIRS} ")

# --------------------------------------------------
# Zoltan libs for plugin
# --------------------------------------------------
if(PV_ZOLTAN_USE_ZOLTAN_1)
  SET(TRILINOS_LIBS zoltan)
else()
  SET(TRILINOS_LIBS zoltan zoltan2)
endif()

# --------------------------------------------------
# Zoltan libs for plugin
# --------------------------------------------------

if(NOT PV_ZOLTAN_USE_ZOLTAN_1)
  SET(ZOLTAN_V2_FILTER
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanV2PartitionFilter.cxx
  )
  SET(ZOLTAN_V2_XML
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanV2PartitionFilter.xml
  )
endif()

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
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanV1PartitionFilter.xml
    ${ZOLTAN_V2_XML}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMeshPartitionFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPartitionOutline.xml

  SERVER_SOURCES
    ${Module_SRCS}
    ${Module_CUDA_C_SRCS}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPKdTree2.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanBasePartitionFilter.cxx

  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkZoltanV1PartitionFilter.cxx
    ${ZOLTAN_V2_FILTER}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMeshPartitionFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkBoundsExtentTranslator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPartitionOutline.cxx

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
# Enable Testing
#--------------------------------------------------
include( CTest )
enable_testing()

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
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src/algorithms"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src/environment"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src/input"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src/models"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src/problems"
  "${PROJECT_SOURCE_DIR}/trilinos/packages/zoltan2/src/util"
)
if(WIN32)
  set(${PROJECT_NAME}_INCLUDES_BUILD_TIME 
    ${PROJECT_SOURCE_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include
    ${${PROJECT_NAME}_INCLUDES_BUILD_TIME}      
  )
endif()


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
