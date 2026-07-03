configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.hxx.in"
  "${CMAKE_CURRENT_BINARY_DIR}/vtkParticlePartitionFilter.h"
  @ONLY)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/vtkMeshPartitionFilter.hxx.in"
  "${CMAKE_CURRENT_BINARY_DIR}/vtkMeshPartitionFilter.h"
  @ONLY)

set(sources
  vtkZoltanBasePartitionFilter.cxx
  vtkZoltanV2PartitionFilter.cxx
  vtkParticlePartitionFilter.cxx
  vtkMeshPartitionFilter.cxx
  vtkPartitionOutline.cxx
  vtkPKdTree2.cxx
  vtkZoltanBasePartitionFilter.txx)

set(headers
  vtkZoltanBasePartitionFilter.h
  vtkZoltanV2PartitionFilter.h
  vtkZoltanVersion.h
  vtkPartitionOutline.h
  vtkPKdTree2.h
  "${CMAKE_CURRENT_BINARY_DIR}/vtkParticlePartitionFilter.h"
  "${CMAKE_CURRENT_BINARY_DIR}/vtkMeshPartitionFilter.h")

vtk_module_add_module(pvZoltanFilters
  SOURCES ${sources}
  HEADERS ${headers})



# System Trilinos depends on Kokkos but may not propagate its include path.
find_path(Kokkos_INCLUDE_DIR
  NAMES Kokkos_Core.hpp
  PATHS "$ENV{KOKKOS_ROOT}/include"
        "/home/biddisco/src/_env/spack/environments/opal/spack-env/view/kokkos/5.1.1/include"
  DOC "Kokkos include directory")

# Trilinos also depends on KokkosKernels; its include path is not exported.
find_path(KokkosKernels_INCLUDE_DIR
  NAMES KokkosSparse_Utils.hpp
  PATHS "$ENV{KOKKOSKERNELS_ROOT}/include"
        "/home/biddisco/src/_env/spack/environments/opal/spack-env/view/kokkos-kernels/5.1.1/include"
  DOC "KokkosKernels include directory")

# ParMETIS/METIS headers are needed by Zoltan2 when ParMETIS is enabled.
find_path(ParMETIS_INCLUDE_DIR
  NAMES parmetis.h
  PATHS "$ENV{PARMETIS_ROOT}/include"
        "/home/biddisco/src/_env/spack/environments/opal/spack-env/view/parmetis/4.0.3/include"
  DOC "ParMETIS include directory")

find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  PATHS "$ENV{METIS_ROOT}/include"
        "/home/biddisco/src/_env/spack/environments/opal/spack-env/view/metis/5.1.0/include"
  DOC "METIS include directory")

# Kokkos with CUDA support needs the CUDA toolkit headers.
find_package(CUDAToolkit REQUIRED)

# System Trilinos does not export Kokkos link libraries either.
find_library(ZOLTAN2_LIBRARY
  NAMES zoltan2
  PATHS "${Trilinos_INSTALL_DIR}/lib"
        "/home/biddisco/opt/spack.git/opt/spack/linux-skylake/trilinos-17.1.1-2yinozmcwuhz7c4q7hwn74fk5vel23l5/lib"
  DOC "Zoltan2 library")

find_library(KOKKOS_CORE_LIBRARY
  NAMES kokkoscore
  PATHS "$ENV{KOKKOS_ROOT}/lib"
        "/home/biddisco/src/_env/spack/environments/opal/spack-env/view/kokkos/5.1.1/lib"
  DOC "Kokkos core library")

find_library(KOKKOS_CONTAINERS_LIBRARY
  NAMES kokkoscontainers
  PATHS "$ENV{KOKKOS_ROOT}/lib"
        "/home/biddisco/src/_env/spack/environments/opal/spack-env/view/kokkos/5.1.1/lib"
  DOC "Kokkos containers library")

vtk_module_include(pvZoltanFilters
  PUBLIC
    "$<BUILD_INTERFACE:${Trilinos_INCLUDE_DIRS}>"
    "$<BUILD_INTERFACE:${Kokkos_INCLUDE_DIR}>"
    "$<BUILD_INTERFACE:${KokkosKernels_INCLUDE_DIR}>"
    "$<BUILD_INTERFACE:${ParMETIS_INCLUDE_DIR}>"
    "$<BUILD_INTERFACE:${METIS_INCLUDE_DIR}>"
    "$<BUILD_INTERFACE:${CUDAToolkit_INCLUDE_DIRS}>")

# Use long long global IDs to match the Trilinos/Tpetra explicit instantiations.
vtk_module_definitions(pvZoltanFilters PUBLIC HAVE_ZOLTAN2_LONG_LONG_INT)

# The base filter source still uses the legacy VTK_USE_MPI guard.
vtk_module_definitions(pvZoltanFilters PUBLIC VTK_USE_MPI)

# Trilinos headers pull in Kokkos; tell the Kokkos compiler wrapper to route
# these translation units through nvcc_wrapper.
vtk_module_definitions(pvZoltanFilters PRIVATE KOKKOS_DEPENDENCE)

# Kokkos/Trilinos use __host__/__device__ lambdas.
vtk_module_compile_options(pvZoltanFilters PRIVATE --extended-lambda)

# Ensure the actual shared library links against Trilinos/Kokkos.
# pvZoltanFilters is an ALIAS; resolve it before calling target_link_libraries.
get_target_property(_pvZoltanFilters_real pvZoltanFilters ALIASED_TARGET)
target_link_libraries(${_pvZoltanFilters_real}
  PUBLIC
    ${ZOLTAN2_LIBRARY}
    ${KOKKOS_CORE_LIBRARY}
    ${KOKKOS_CONTAINERS_LIBRARY})

paraview_add_server_manager_xmls(
  XMLS
    vtkZoltanBasePartitionFilter.xml
    vtkZoltanV2PartitionFilter.xml
    vtkParticlePartitionFilter.xml
    vtkMeshPartitionFilter.xml
    vtkPartitionOutline.xml)
