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
# a macro to add one test
#--------------------------------------------------
MACRO(add_local_test test_name) 
  ADD_EXECUTABLE(${test_name} ${test_name}.cxx)
  TARGET_LINK_LIBRARIES(${test_name}
    pv_meshless 
    )
  ADD_TEST(${test_name} ${PROJECT_BINARY_DIR}/../bin/${test_name})
ENDMACRO(add_local_test)

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
# Simple non MPI Filter/Kernel Tests 
#--------------------------------------------------
SET(TEST_LIST 
)

FOREACH(test_name ${TEST_LIST})
  add_local_test(${test_name})
ENDFOREACH(test_name)

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
    pv_zoltan 
    ${Trilinos_LIBRARIES}
    )
  
  greater_equal(${MPIEXEC_MAX_NUMPROCS} 4 processors)
  SET(test_name "TestParticlePartitionFilter-P4")
  IF (processors)
    ADD_TEST(
      NAME ${test_name}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
        $<TARGET_FILE:TestParticlePartitionFilter> 
        -testName ${test_name}
        -generateParticles 1000
        -V ${PROJECT_SOURCE_DIR}/Testing/baseline/ParticlePartition.png
        -doRender 1
        -windowSize "400 400"
    )
  ENDIF (processors)
 

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
