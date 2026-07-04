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
    NAME TestParticlePartitionWeightCount-P4-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionWeightCount>
      -testName TestParticlePartitionWeightCount-P4
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 5000
      -particleGenerator 1
      -useWeights 1
      -scalarName ProcessId
      -scalarMode 0
      -scalarRange "0 3"
      -actorShift 0.5
      -doRender 0
      -windowSize "800 800"
  )

  #------------------------------------------------
  # Particle partition tests with ghost layers
  #------------------------------------------------
  add_test(
    NAME TestParticlePartitionFilterGhost-P4-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter>
      -testName TestParticlePartitionFilterGhost-P4
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 1000
      -particleGenerator 0
      -ghostOverlap 100
      -scalarName vtkGhostType
      -scalarMode 0
      -scalarRange "0 1"
      -actorShift 0.5
      -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/ParticlePartitionGhosts-${_test_version}.png
      -doRender 1
      -windowSize "800 800"
  )

  #------------------------------------------------
  # Particle partition tests with weights (rendered)
  #------------------------------------------------
  add_test(
    NAME TestParticlePartitionWeights-P4-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter>
      -testName TestParticlePartitionWeights-P4
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 5000
      -particleGenerator 1
      -ghostOverlap 50
      -scalarName Weights
      -scalarMode 0
      -scalarRange "0 1"
      -actorShift 0.5
      -useWeights 1
      -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/ParticlePartitionWeights-${_test_version}.png
      -doRender 1
      -windowSize "800 800"
  )

  #------------------------------------------------
  # Mesh partition tests: boundary modes
  #------------------------------------------------
  set(boundary_list "first;most;centroid")
  foreach(dmode soma sphere)
    # Sphere variants need a slightly higher threshold due to driver/AA differences.
    if(dmode STREQUAL "sphere")
      set(threshold_ 0)
    else()
      set(threshold_ 0)
    endif()
    foreach(mode 0 1 2)
      list(GET boundary_list ${mode} BMODE)
      add_test(
        NAME TestMeshPartitionFilterBoundary-${BMODE}-${dmode}-${_test_version}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
          $<TARGET_FILE:TestMeshPartitionFilter>
          -testName TestMeshPartitionFilterBoundary-${BMODE}-${dmode}
          -T "${PLUGIN_TEST_DIR}"
          -F ${dmode}.vtp
          -D ${CMAKE_CURRENT_SOURCE_DIR}/data
          -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartition-${BMODE}-${dmode}-${_test_version}.png
          -doRender 1
          -edges 1
          -ghostMode 0
          -boundaryMode ${mode}
          -scalarMode 0
          -scalarName ProcessId
          -scalarRange "0 3"
          -actorShift 0.0
          -imageThreshold ${threshold_}
          -windowSize "1024 800"
      )
    endforeach()
  endforeach()

  #------------------------------------------------
  # Mesh partition tests: cell data scalars
  #------------------------------------------------
  set(scalar_vars "RTNeuron Opacity" "Elevation")
  set(index 0)
  foreach(dmode soma sphere)
    list(GET scalar_vars ${index} SCALAR)
    add_test(
      NAME TestMeshPartitionFilterCellData-${dmode}-${_test_version}
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
        $<TARGET_FILE:TestMeshPartitionFilter>
        -testName TestMeshPartitionFilterCellData-${dmode}
        -T "${PLUGIN_TEST_DIR}"
        -F ${dmode}.vtp
        -D ${CMAKE_CURRENT_SOURCE_DIR}/data
        -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartition-celldata-${dmode}-${_test_version}.png
        -doRender 1
        -ghostMode 0
        -scalarMode 1
        -scalarName "${SCALAR}"
        -scalarRange "0 1"
        -actorShift 0.0
        -windowSize "1024 800"
    )
    math(EXPR index "${index} + 1")
  endforeach()

  #------------------------------------------------
  # Mesh partition tests: ghost modes
  #------------------------------------------------
  set(ghostmode_list "none;boundary;bbox")
  foreach(dmode soma sphere)
    if(dmode STREQUAL "sphere")
      set(threshold_ 0)
    else()
      set(threshold_ 0)
    endif()
    set(index 0)
    foreach(ghostmode 0 1 2)
      list(GET ghostmode_list ${ghostmode} GMODE)
      set(overlap_ "")
      set(actorshift_ -actorShift 0.25)
      if(ghostmode STREQUAL "2")
        if(dmode STREQUAL "sphere")
          set(overlap_ -ghostOverlap 0.09)
          set(actorshift_ -actorShift 0.6)
        else()
          set(overlap_ -ghostOverlap 2.5)
          set(actorshift_ -actorShift 0.5)
        endif()
      endif()

      add_test(
        NAME TestMeshPartitionFilterGhost-${GMODE}-${dmode}-${_test_version}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
          $<TARGET_FILE:TestMeshPartitionFilter>
          -testName TestMeshPartitionFilterGhost-${GMODE}-${dmode}
          -T "${PLUGIN_TEST_DIR}"
          -F ${dmode}.vtp
          -D ${CMAKE_CURRENT_SOURCE_DIR}/data
          -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartition-ghost-${GMODE}-${dmode}-${_test_version}.png
          -doRender 1
          -edges 1
          -boundaryMode 2
          -ghostMode ${index}
          ${overlap_}
          -scalarMode 1
          -scalarName vtkGhostRanks
          -scalarRange "0 4"
          ${actorshift_}
          -imageThreshold ${threshold_}
          -windowSize "1024 800"
      )
      math(EXPR index "${index} + 1")
    endforeach()
  endforeach()

  #------------------------------------------------
  # Mesh partition tests: unstructured grid
  #------------------------------------------------
  set(threshold_ 0)
  set(ghostmode_list "none;boundary;bbox")
  set(index 0)
  foreach(ghostmode 0 1 2)
    list(GET ghostmode_list ${ghostmode} GMODE)
    set(overlap_ "")
    if(ghostmode STREQUAL "0")
      set(scalars_ -scalarMode 0 -scalarName RTData -scalarRange "0 1")
    elseif(ghostmode STREQUAL "1")
      set(scalars_ -scalarMode 1 -scalarName vtkGhostRanks -scalarRange "0 4")
    elseif(ghostmode STREQUAL "2")
      set(overlap_ -ghostOverlap 0.2)
      set(scalars_ -scalarMode 1 -scalarName vtkGhostRanks -scalarRange "0 4")
    endif()

    add_test(
      NAME TestMeshPartitionFilterUnstructured-${GMODE}-${_test_version}
      COMMAND
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
        $<TARGET_FILE:TestMeshPartitionFilter>
        -testName TestMeshPartitionFilterUnstructured-${GMODE}
        -T "${PLUGIN_TEST_DIR}"
        -F wavelet-unstructured.vtu
        -D ${CMAKE_CURRENT_SOURCE_DIR}/data
        -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartition-unstructured-${GMODE}-${_test_version}.png
        -doRender 1
        -edges 1
        -boundaryMode 2
        -ghostMode ${index}
        ${overlap_}
        ${scalars_}
        -actorShift 1.0
        -imageThreshold ${threshold_}
        -windowSize "1024 800"
    )
    math(EXPR index "${index} + 1")
  endforeach()

  #------------------------------------------------
  # Mesh partition tests: point-data scalar transfer
  #------------------------------------------------
  add_test(
    NAME TestMeshPartitionFilterScalars-P4-${_test_version}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestMeshPartitionFilter>
      -testName TestMeshPartitionFilterScalars-P4
      -T "${PLUGIN_TEST_DIR}"
      -F soma-4-partitions.pvtp
      -D ${CMAKE_CURRENT_SOURCE_DIR}/data
      -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartitionScalars-${_test_version}.png
      -boundaryMode 1
      -doRender 1
      -scalarMode 0
      -scalarName "RTNeuron Opacity"
      -scalarRange "0 1"
      -actorShift 0.25
      -windowSize "1024 800"
  )

  #------------------------------------------------
  # Mesh partition tests: mixed polydata
  #------------------------------------------------
  # Disabled for now: mixed cell types (verts/lines/polys) are not handled
  # correctly by the current mesh partitioner, producing an almost blank image.
  # add_test(
  #   NAME TestMeshPartitionFilterMixedPolys-P4-${_test_version}
  #   COMMAND
  #     ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
  #     $<TARGET_FILE:TestMeshPartitionFilter>
  #     -testName TestMeshPartitionFilterMixedPolys-P4
  #     -T "${PLUGIN_TEST_DIR}"
  #     -F polydata-mixed.pvtp
  #     -D ${CMAKE_CURRENT_SOURCE_DIR}/data
  #     -V ${CMAKE_CURRENT_SOURCE_DIR}/baseline/MeshPartition-mixed-polydata.png
  #     -boundaryMode 1
  #     -doRender 1
  #     -cameraPosition "50 140 50"
  #     -cameraFocus    "12 0 12"
  #     -cameraViewUp   "-0.6 0.4 -0.6"
  #     -scalarMode 0
  #     -scalarName "Elevation"
  #     -scalarRange "0 1"
  #     -actorShift 0.25
  #     -windowSize "1024 800"
  # )
endif ()
