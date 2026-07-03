if (NOT BUILD_TESTING)
  return()
endif ()

set(PLUGIN_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary")
file(MAKE_DIRECTORY "${PLUGIN_TEST_DIR}")

find_package(MPI REQUIRED)
find_package(CUDAToolkit REQUIRED)

set(_test_srcs
  TestParticlePartitionFilter.cxx
  TestParticlePartitionWeightCount.cxx
  TestMeshPartitionFilter.cxx
  TestUtils.cxx
)

set(_test_executables
  TestParticlePartitionFilter
  TestParticlePartitionWeightCount
  TestMeshPartitionFilter
)

foreach (_exe IN LISTS _test_executables)
  add_executable(${_exe} ${_exe}.cxx TestUtils.cxx)
  target_link_libraries(${_exe}
    PRIVATE
      VTK::TestingRendering
      pvZoltanFilters
      Zoltan2::all_libs
      MPI::MPI_CXX
      CUDA::cudart
  )
  vtk_module_autoinit(
    TARGETS ${_exe}
    MODULES VTK::TestingRendering pvZoltanFilters
  )

  # Test utilities call Kokkos directly; route these TUs through the Kokkos
  # compiler wrapper just like the filter library.
  target_compile_definitions(${_exe}
    PRIVATE
      KOKKOS_DEPENDENCE)
  target_compile_options(${_exe}
    PRIVATE
      --extended-lambda)
endforeach ()

if (MPIEXEC_MAX_NUMPROCS LESS 4)
  message(WARNING "pv-zoltan tests are designed to run with at least 4 MPI ranks")
endif ()

if (MPIEXEC_MAX_NUMPROCS GREATER_EQUAL 4)
  set(_test_version "v2")

  add_test(
    NAME TestParticlePartitionFilter-P4-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter>
      -testName TestParticlePartitionFilter-P4
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 1000
      -particleGenerator 0
      -scalarName ProcessId
      -scalarMode 0
      -scalarRange "0 3"
      -actorShift 0.5
      -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/ParticlePartition-${_test_version}.png
      -doRender 1
      -windowSize "800 800"
  )

  add_test(
    NAME TestMeshPartitionFilterBoundary-first-soma-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestMeshPartitionFilter>
      -testName TestMeshPartitionFilterBoundary-first-soma
      -T "${PLUGIN_TEST_DIR}"
      -F soma.vtp
      -D ${CMAKE_CURRENT_SOURCE_DIR}/data
      -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartition-first-soma-${_test_version}.png
      -doRender 1
      -edges 1
      -ghostMode 0
      -boundaryMode 0
      -scalarMode 0
      -scalarName ProcessId
      -scalarRange "0 3"
      -actorShift 0.0
      -windowSize "1024 800"
  )

  add_test(
    NAME TestParticlePartitionWeightCount-P4-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionWeightCount>
      -testName TestParticlePartitionWeightCount-P4
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 300
      -particleGenerator 1
      -useWeights 1
      -scalarName ProcessId
      -scalarMode 0
      -scalarRange "0 3"
      -actorShift 0.5
      -doRender 0
      -windowSize "800 800"
  )
endif ()
