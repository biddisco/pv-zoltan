project("pv_zoltan")

#--------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
endif(COMMAND cmake_policy)

add_definitions(-DCMAKE_IS_BROKEN)


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

#--------------------------------------------------
# Find cuda for optional Piston classes
#--------------------------------------------------
option(PV_ZOLTAN_USE_CUDA "Use CUDA rendering" OFF)

if(PV_ZOLTAN_USE_CUDA)

  #--------------------------------------------------
  # setup module definitions we'll need for piston
  # compiler flags/cuda etc
  #--------------------------------------------------
  include(vtkModuleAPI)
  include(vtkModuleMacros)
  include(ParaViewMacros)

  #--------------------------------------------------
  # Find CUDA and Thrust
  #--------------------------------------------------
  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/piston")
  find_package(CUDA REQUIRED)
  find_package(CUDAThrust REQUIRED)

  #--------------------------------------------------
  # Find piston header library
  #--------------------------------------------------
  find_path(PISTON_INCLUDE
    PistonCodingStyle.xml
    "/Source/PISTON/github/PISTON"
    DOC "PISTON library's source directory")

  if(PISTON_INCLUDE)
    add_definitions(-DPV_ZOLTAN_USE_PISTON)

    INCLUDE_DIRECTORIES(
      ${CMAKE_CURRENT_SOURCE_DIR}/piston
    )
  endif(PISTON_INCLUDE)

endif(PV_ZOLTAN_USE_CUDA)

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
# Trilinos support 
# --------------------------------------------------
SET(TRILINOS_ENABLED 0)

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
#  ADD_DEFINITIONS(-DREALLOC_BUG)
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

#--------------------------------------------------
# Define PISTON sources
#--------------------------------------------------
if(PV_ZOLTAN_USE_CUDA)
  set(Module_SRCS
    piston/vtkDataSetToPiston.cxx
    piston/vtkPistonAlgorithm.cxx
    piston/vtkPistonDataObject.cxx
    piston/vtkPistonPolygonsPainter.cxx
  )

  set(Module_CUDA_SRCS
    piston/vtkPistonConverters.cu
    piston/vtkPistonMapper.cu
  )

  #--------------------------------------------------
  # include path for Thrust
  #--------------------------------------------------
  include_directories(${CUDATHRUST_INCLUDE_DIRS})
  message("include_directories(${CUDATHRUST_INCLUDE_DIRS})")

  #--------------------------------------------------
  # Define CUDA sources
  #--------------------------------------------------
  if(BUILD_SHARED_LIBS)
    set(GPGPU_BUILD_TYPE SHARED)
  endif()

  # set this, otherwise C_SOURCES is empty in VS
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF) 
  # generate the list of dependencies based on our module definition file
  vtk_module_impl()
  # remove a bunch of cruft that VTK adds in here about modules autoinit
  # first get the old setting and store it
  get_directory_property(TEMP_INCLUDE_DIRECTORIES INCLUDE_DIRECTORIES)
  get_directory_property(TEMP_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)
  # clear compile defs 
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "")
  # shorten the include dirs so nvcc doesn't crumble under the entire VTK/ParaView/Qt path
  set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${${vtk-module}_DEPENDS_INCLUDE_DIRS};${PISTON_INCLUDE}")
  # call the cuda wrapper
  cuda_compile(Module_CUDA_C_SRCS ${Module_CUDA_SRCS} ${GPGPU_BUILD_TYPE})
  # put the include dirs back to what they were
  set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${TEMP_INCLUDE_DIRECTORIES}")
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${TEMP_COMPILE_DEFINITIONS}")

  set_source_files_properties(
    ${Module_CUDA_C_SRCS}
    WRAP_EXCLUDE
  )
endif(PV_ZOLTAN_USE_CUDA)

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

    ${XML_PATH}/vtkDepthSortRepresentation.xml

  SERVER_MANAGER_SOURCES
    ${TRILINOS_CXX}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkBoundsExtentTranslator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPartitionOutline.cxx

    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortDefaultPainter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortPainter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortPolyData2.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortPolygonsPainter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTwoScalarsToColorsPainter.cxx

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

if(CMAKE_COMPILER_IS_GNUCXX)
  message("Adding C++Ox flag for Depth Sort")
  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/vtkDepthSortPolyData2.cxx PROPERTIES COMPILE_FLAGS "-std=c++0x")
endif(CMAKE_COMPILER_IS_GNUCXX)

#--------------------------------------------------
# setup plugin linking
#--------------------------------------------------
TARGET_LINK_LIBRARIES(${PLUGIN_NAME} 
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
    ${CUDA_LIBRARIES}
)

#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(testing)
ENDIF (BUILD_TESTING)  


