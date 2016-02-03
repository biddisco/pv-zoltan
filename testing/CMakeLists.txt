#--------------------------------------------------
# Tests here are added using the 
# ADD_TEST(NAME ... COMMAND ...)
# form and will only run on windows when ctest is run using
# ctest -C Release / Debug / etc
#--------------------------------------------------

#--------------------------------------------------
# Machines running OSMesa instead of OpenGL might need a larger threshold for image test comparisons
#--------------------------------------------------
set(IMAGE_THRESHOLD "0")

#--------------------------------------------------
# convenience macro for a>=b for testing processor count
#--------------------------------------------------
MACRO(greater_equal a b boolresult) 
  IF ((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
    SET(${boolresult} 1)
  ELSE((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
    SET(${boolresult} 0)
  ENDIF ((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
ENDMACRO(greater_equal)

#--------------------------------------------------
# Slurm job manager require "-rmk" "slurm" arguments
# so split MPIEXEC_NUMPROC_FLAG into a list of args
# to stop cmake from putting quotes around all of it
#--------------------------------------------------
string(REPLACE " " ";" MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG})

#--------------------------------------------------
# Trilinos
#--------------------------------------------------
ADD_EXECUTABLE(TestParticlePartitionFilter 
  TestParticlePartitionFilter.cxx
  TestUtils.cxx
)
TARGET_LINK_LIBRARIES(TestParticlePartitionFilter
  vtkTestingRendering
  pv_zoltan 
  ${Trilinos_LIBRARIES}
)

ADD_EXECUTABLE(TestParticlePartitionWeightCount
  TestParticlePartitionWeightCount.cxx
  TestUtils.cxx
)
TARGET_LINK_LIBRARIES(TestParticlePartitionWeightCount
  vtkTestingRendering
  pv_zoltan 
  ${Trilinos_LIBRARIES}
)

ADD_EXECUTABLE(TestMeshPartitionFilter 
  TestMeshPartitionFilter.cxx
  TestUtils.cxx
)
TARGET_LINK_LIBRARIES(TestMeshPartitionFilter
  vtkTestingRendering
  pv_zoltan 
  ${Trilinos_LIBRARIES}
)
  
#ADD_EXECUTABLE(TestPartitionCompositing
#  TestPartitionCompositing.cxx
#  TestUtils.cxx
#)
#TARGET_LINK_LIBRARIES(TestPartitionCompositing
#  pv_zoltan 
#  ${Trilinos_LIBRARIES}
#)

#--------------------------------------------------
# Tests will only run if we have 4+ processors
#--------------------------------------------------

greater_equal(${MPIEXEC_MAX_NUMPROCS} 4 processors)
if (NOT processors)
    message(WARNING "pv-zoltan tests require at least 4 mpi ranks in MPIEXEC_MAX_NUMPROCS") 
endif()

IF (processors)

  #------------------------------------------------
  # Particle partition tests
  #------------------------------------------------
  SET(test_name "TestParticlePartitionFilter-P4")
  ADD_TEST(
    NAME ${test_name}
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter> 
      -testName ${test_name}
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 1000
      -particleGenerator 0
      -scalarName ProcessId
      -scalarMode 0
      -scalarRange "0 3"
      -actorShift 0.5
      -V ${PROJECT_SOURCE_DIR}/testing/baseline/ParticlePartition.png
      -doRender 1
      -windowSize "800 800"
  )
  
  SET(test_name "TestParticlePartitionFilterGhost-P4")
  ADD_TEST(
    NAME ${test_name}
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter> 
      -testName ${test_name}
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 1000
      -particleGenerator 0
      -ghostOverlap 100
      -scalarName vtkGhostLevels
      -scalarMode 0
      -scalarRange "0 1"
      -actorShift 0.5
      -V ${PROJECT_SOURCE_DIR}/testing/baseline/ParticlePartitionGhosts.png
      -doRender 1
      -windowSize "800 800"
  )
  
  SET(test_name "TestParticlePartitionWeights-P4")
  ADD_TEST(
    NAME ${test_name}
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter> 
      -testName ${test_name}
      -T "${PLUGIN_TEST_DIR}"
      -generateParticles 5000
      -particleGenerator 1
      -ghostOverlap 50
      -scalarName Weights
      -scalarMode 0
      -scalarRange "0 1"
      -actorShift 0.5
      -useWeights 1
      -V ${PROJECT_SOURCE_DIR}/testing/baseline/ParticlePartitionWeights.png
      -doRender 1
      -windowSize "800 800"
  )
  
  SET(test_name "TestParticlePartitionWeightCount-P4")
  ADD_TEST(
    NAME ${test_name}
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionWeightCount> 
      -testName ${test_name}
      -generateParticles 5000
      -particleGenerator 1
      -useWeights 1
  )
  
  #------------------------------------------------
  # Mesh partition tests
  #------------------------------------------------
  SET(test_name "TestMeshPartitionFilterBoundary-P4")
  set(boundary_list "first;most;centroid")
  foreach(dmode soma sphere)
      foreach(mode 0 1 2)
          list(GET "boundary_list" "${mode}" BMODE)
          ADD_TEST(
            NAME ${test_name}-${BMODE}-${dmode}
            COMMAND 
              ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
              $<TARGET_FILE:TestMeshPartitionFilter> 
              -testName ${test_name}-${BMODE}-${dmode}
              -T "${PLUGIN_TEST_DIR}"
              -F ${dmode}-4-partitions.pvtp
              -D ${PROJECT_SOURCE_DIR}/testing/data
              -V ${PROJECT_SOURCE_DIR}/testing/baseline/MeshPartition-${BMODE}-${dmode}.png
              -doRender 1
              -ghostMode 0
              -boundaryMode ${mode}
              -scalarMode 0
              -scalarName ProcessId
              -scalarRange "0 3"
              -actorShift 0.0
              -windowSize "1024 800"
          )
      endforeach()
  endforeach()

  SET(test_name "TestMeshPartitionFilterCellData-P4")
  set(scalar_vars "RTNeuron Opacity" "Elevation")
  set(index 0)
  foreach(dmode soma sphere)
    list(GET "scalar_vars" "${index}" SCALAR)
    ADD_TEST(
      NAME ${test_name}-${dmode}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
        $<TARGET_FILE:TestMeshPartitionFilter> 
        -testName ${test_name}-${dmode}
        -T "${PLUGIN_TEST_DIR}"
        -F ${dmode}-4-partitions.pvtp
        -D ${PROJECT_SOURCE_DIR}/testing/data
        -V ${PROJECT_SOURCE_DIR}/testing/baseline/MeshPartition-celldata-${dmode}.png
        -doRender 1
        -ghostMode 0
        -scalarMode 1
        -scalarName ${SCALAR}
        -scalarRange "0 1"
        -actorShift 0.0
        -windowSize "1024 800"
    )
    math(EXPR index "${index} + 1")
  endforeach()

  set(test_name "TestMeshPartitionFilterGhost-P4")
  set(ghostmode_list "none;boundary;bbox;neighbour")
  foreach(dmode soma sphere)
      set(index 0)
      foreach(ghostmode 0 1 2 3)
          list(GET "ghostmode_list" "${ghostmode}" GMODE)
          set(overlap_ "")
          #
          if ("${dmode}" STREQUAL "sphere" )
            set(actorshift_ -actorShift 0.25)
          else()
            set(actorshift_ -actorShift 0.25)
          endif()

          if("${ghostmode}" STREQUAL "2")
            if("${dmode}" STREQUAL "sphere")
              set(overlap_ -ghostOverlap 0.09)
              set(actorshift_ -actorShift 0.6)
            else()
              set(overlap_ -ghostOverlap 2.5)
              set(actorshift_ -actorShift 0.5)
            endif()
          endif()

          ADD_TEST(
            NAME ${test_name}-${GMODE}-${dmode}
            COMMAND 
              ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
              $<TARGET_FILE:TestMeshPartitionFilter> 
              -testName ${test_name}-${GMODE}-${dmode}
              -T "${PLUGIN_TEST_DIR}"
              -F ${dmode}-4-partitions.pvtp
              -D ${PROJECT_SOURCE_DIR}/testing/data
              -V ${PROJECT_SOURCE_DIR}/testing/baseline/MeshPartition-ghost-${GMODE}-${dmode}.png
              -doRender 1
              -edges 1
              -boundaryMode 2
              -ghostMode ${index}
              ${overlap_}
              -scalarMode 1
              -scalarName vtkGhostRanks
              -scalarRange "0 4"
              ${actorshift_}
              -windowSize "1024 800"
          )
          math(EXPR index "${index} + 1")
      endforeach()
  endforeach()



  SET(test_name "TestMeshPartitionFilterUnstructured-P4")
  ADD_TEST(
    NAME ${test_name}
    COMMAND
    ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
    $<TARGET_FILE:TestMeshPartitionFilter>
    -testName ${test_name}
    -T "${PLUGIN_TEST_DIR}"
    -F wavelet-unstructured.vtu
    -D ${PROJECT_SOURCE_DIR}/testing/data
    -V ${PROJECT_SOURCE_DIR}/testing/baseline/MeshPartitionScalars.png
    -boundaryMode 1
    -doRender 1
    -scalarMode 0
    -scalarName "RTData"
    -scalarRange "0 1"
    -actorShift 0.25
    -windowSize "1024 800"
  )

  SET(test_name "TestMeshPartitionFilterScalars-P4")
  ADD_TEST(
    NAME ${test_name}
    COMMAND
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestMeshPartitionFilter>
      -testName ${test_name}
      -T "${PLUGIN_TEST_DIR}"
      -F soma-4-partitions.pvtp
      -D ${PROJECT_SOURCE_DIR}/testing/data
      -V ${PROJECT_SOURCE_DIR}/testing/baseline/MeshPartitionScalars.png
      -boundaryMode 1
      -doRender 1
      -scalarMode 0
      -scalarName "RTNeuron Opacity"
      -scalarRange "0 1"
      -actorShift 0.25
      -windowSize "1024 800"
  )
ENDIF (processors)
 
#SET(test_name "TestPartitionCompositing-P4")
#IF (processors)
#  ADD_TEST(
#    NAME ${test_name}
#    COMMAND 
#      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
#      $<TARGET_FILE:TestPartitionCompositing> 
#      -testName ${test_name}
#      -T "${PLUGIN_TEST_DIR}"
#      -F neuron-test.vtp
#      -D ${PROJECT_SOURCE_DIR}/testing/data
#      -V ${PROJECT_SOURCE_DIR}/testing/baseline/MeshPartition.png
#      -doRender 1
#      -cameraPosition "-100.613214826891 1393.997094168 353.866251068379"
#      -cameraFocus    "-100.613214826891 1393.997094168 -1565.74007414479"
#      -cameraViewUp   "0 1 0"
#      -scalar "vtkProcessId"
#      -windowSize "1024 800"
#  )
#ENDIF (processors)

#--------------------------------------------------
#--------------------------------------------------
# Benchmark SC12 : Take one dataset with bounding boxes
# saved and run partitioning filter on different CPU counts
# with/without boxes. Time the partitioning.
#--------------------------------------------------
#--------------------------------------------------
LIST(GET SC12_DATA_BOXN 0 SC12_TEST)
IF (EXISTS ${SC12_DATA_PATH}/${SC12_DATA_FILE}${SC12_TEST}.h5part)
  FOREACH(index 0 1)
    SET(ncount 0)
    FOREACH(N ${SC12_DATA_BOXN})
      LIST(GET SC12_DATA_DX_N ${ncount} DELTAX)
      LIST(GET SC12_DATA_DX_4 ${ncount} GST)
      SET(TestParams 
        -D "${SC12_DATA_PATH}"
        -F "${SC12_DATA_FILE}${N}.h5part"
        -T "${PLUGIN_TEST_DIR}"
        -Xarray "X"
        -Yarray "Y"
        -Zarray "Z"
        -particlesize ${DELTAX}
        -ghost_region ${GST}
        -ignorePartitions "${index}"
        -randomizeExtents "${index}"
      )
      # remove 000 from 00032 for Nprocs etc
      STRING(REGEX REPLACE "(^0+)" "" N "${N}")
      #
      SET(M ${N})
      while (NOT ${M} LESS ${SC12_DATA_NMIN})
        IF (NOT ${M} GREATER ${MPIEXEC_MAX_NUMPROCS})
          SET(test_name "benchmark-PartitionBoxes-B${N}-P${M}-I${index}")
          MATH(EXPR MPI_DEPTH "(${N}/${CORES_PER_NODE})/${M}")
          IF (${MPI_DEPTH} GREATER 1)
            SET (MPI_DEPTH_FLAGS -d ${MPI_DEPTH})
          ELSE (${MPI_DEPTH} GREATER 1)
            SET (MPI_DEPTH_FLAGS )
          ENDIF (${MPI_DEPTH} GREATER 1)
          IF (${MPI_DEPTH} GREATER ${CORES_PER_NODE})
            SET (MPI_DEPTH_FLAGS -d ${CORES_PER_NODE})
          ENDIF (${MPI_DEPTH} GREATER ${CORES_PER_NODE})
          ADD_TEST(
            NAME ${test_name}
            COMMAND 
              ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${M} ${MPI_DEPTH_FLAGS}
              $<TARGET_FILE:TestParticlePartitionFilter> 
              -testName ${test_name}
              ${TestParams}
          )
        ENDIF (NOT ${M} GREATER ${MPIEXEC_MAX_NUMPROCS})
        #
        MATH(EXPR M "${M} / 2")
      endwhile (NOT ${M} LESS ${SC12_DATA_NMIN})
      MATH(EXPR ncount "${ncount} + 1")
    ENDFOREACH(N)
  ENDFOREACH(index)
ENDIF (EXISTS ${SC12_DATA_PATH}/${SC12_DATA_FILE}${SC12_TEST}.h5part)
