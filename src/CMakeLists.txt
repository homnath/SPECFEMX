# REVISION:
#   WE, Jan 25th 2023; 
#   HNG, Mar 21,2017; HNG, Sep 06,2011;

# The following may be useful for reference
# Libraries: 
#   sem    - core library with SEM source codes, (shared by all exec)
#   serial - library with codes for serial operation  
#   parsem - library with codes for parallel operation  
#   both   - library with codes from 'srcboth' 

# Variables: 
#   srcsem     - array of source file names for SEM library 
#   srcserial  - array of source files for serial library
#   srcparsem  - array of source file names for Parallel (parsem) library 
#   srcboth    - array of source files required by both serial and parallel builds

# Dependencies: 
#   serial depends on SEM 
#   parsem depends on SEM 
#   If in parallel build, then 'both' library depends on sem and parsem 
#   If in serial   build, then 'both' library depends on sem and serial 

project(specfemx 
        VERSION 
        1.0
        )

# Source files for sem library
SET(srcsem
  global.f90
  local.f90
  shared.f90
  output_to_user.F90
  package_version.f90
  string_library.f90
  math_library.f90
  shape_library.f90
  element.f90
  degrees_of_freedom.f90
  elastic.f90
  infinite_layer.f90
  recreate_fault.f90
  read_input.f90
  gll_library.f90
  integration.f90
  infinite_element.f90
  utmgeo.f90
  coordinate.f90
  write_ensight.f90
  fault.f90
  weakform.f90
  prepare_gravity.f90
  apply_bc.f90
  preprocess.f90
  matrix_vector.f90
  gindex.f90
  apply_traction.f90
  apply_mtraction.f90
  source_function.f90
  map_location.f90
  free_surface.f90
  compute_pressure.f90
  viscoelastic.f90
  plastic_library.f90
  okada_solution.f90
  postprocess.f90)
  
# Name of library is 'sem' and contains all of the srcsem files? 
ADD_LIBRARY(sem ${srcsem})


# Build parsem library
IF(ISPARALLEL)
  # Variable name is srcparsem 
  SET(srcparsem
     set_precision_mpi.f90
     mpi_library.f90
     math_library_mpi.f90
     ghost_library_mpi.f90
     prepare_sparse.f90
     parsolver.f90)
     # If using PETSC then add relevant f90 to the list srcparsem
  
  IF(USE_PETSC)
    IF(ISCOMPLEX)
      LIST(APPEND srcparsem
        parsolver_petsc_complex.F90)
    ELSE(ISCOMPLEX)
      LIST(APPEND srcparsem
        parsolver_petsc.F90)
    ENDIF(ISCOMPLEX)
  ENDIF(USE_PETSC)


  # Add a library called parsem with the files held in srcparsem 
  ADD_LIBRARY(parsem ${srcparsem})

  # Add dependencies - Parsem library depends on sem library 
  ADD_DEPENDENCIES(parsem sem)

  # If using PETSC then we need to link the directory in which the petsc libraries are, for PARSEM 
  IF(USE_PETSC)
    TARGET_INCLUDE_DIRECTORIES(parsem PUBLIC "${PETSC_DIR}/include" "${PETSC_ARCH}/include")
    TARGET_LINK_LIBRARIES(parsem ${PETSC_ARCH}/lib/libpetsc.so) # Links this library to the executable parsem 
  ENDIF(USE_PETSC)

# Build serial library
ELSE(ISPARALLEL)
  SET(srcserial
     serial_library.f90
     math_library_serial.f90
     prepare_sparse_serial.f90
     solver.f90
     solver_petsc.f90)
  ADD_LIBRARY(serial ${srcserial})
  ADD_DEPENDENCIES(serial sem)
ENDIF(ISPARALLEL)

SET(srcboth
   count_elements.F90
   relaxation_time.F90
   prestress.F90
   prepare_solver.F90
   mesh_spec.F90
   model.F90
   nondimension.F90
   cmtsolution.F90
   earthquake_load.F90
   electrical_load.F90
   station.F90
   benchmark.F90
   save_mesh.F90
   user_input.F90
   cleanup.F90
   print_w_MPI.F90
   determine_solver.F90
   save_variables.F90
   nonlinearloop.F90
   sea_level.F90
   ice.F90
   time_loop.F90
   bcs_and_dof.F90
   other_forces.F90
   initialise_arrays.F90
   )

# Compiler preprocessor directives
# Set the property of the srcboth variable so that it uses MPI??? 
IF(ISPARALLEL)
  SET_PROPERTY(SOURCE ${srcboth} PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE ${srcboth} PROPERTY INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
ENDIF(ISPARALLEL)


ADD_LIBRARY(both ${srcboth})
IF(ISPARALLEL)
  ADD_DEPENDENCIES(both sem parsem)
ELSE(ISPARALLEL)
  ADD_DEPENDENCIES(both sem serial)
ENDIF(ISPARALLEL)

#IF(ISCOMPLEX)
#ELSE(ISCOMPLEX)
#ENDIF(ISCOMPLEX)


# Build specfem library
#SET(srcspecfem "")
#message("${srcspecfem}")
#IF(BUILD_EXCAVATION)
#  LIST(APPEND srcspecfem semexcav3d.f90)
#ENDIF(BUILD_EXCAVATION)
#IF(BUILD_SLOPE)
#  LIST(APPEND srcspecfem semslope3d.f90)
#ENDIF(BUILD_SLOPE)


# ------------------------ SPECFEM ------------------------

SET(srcspecfem
   main.F90)
# Compiler preprocessor directives
#ADD_DEFINITIONS(CMAKE_Fortran_FLAGS "-DUSE_MPI")
IF(ISPARALLEL)
  #SET_PROPERTY(SOURCE ${srcspecfem} PROPERTY COMPILE_DEFINITIONS USE_MPI)
  SET_PROPERTY(SOURCE ${srcspecfem} APPEND PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE ${srcspecfem} PROPERTY INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
ENDIF(ISPARALLEL)

IF(USE_PETSC)
  SET_PROPERTY(SOURCE ${srcspecfem} APPEND PROPERTY COMPILE_DEFINITIONS USE_PETSC)
  IF(ISCOMPLEX)
    SET_PROPERTY(SOURCE ${srcspecfem} APPEND PROPERTY COMPILE_DEFINITIONS USE_COMPLEX)
  ENDIF(ISCOMPLEX)
ENDIF(USE_PETSC)
ADD_LIBRARY(specfem ${srcspecfem})


# Following PETSC link is necessary due to the PETSC module called in the
# main.F90
IF(USE_PETSC)
  TARGET_INCLUDE_DIRECTORIES(specfem PUBLIC "${PETSC_DIR}/include" "${PETSC_ARCH}/include")
  TARGET_LINK_LIBRARIES(specfem ${PETSC_ARCH}/lib/libpetsc.so)
ENDIF(USE_PETSC)
IF(ISPARALLEL)
  ADD_DEPENDENCIES(specfem sem parsem both)
ELSE(ISPARALLEL)
  ADD_DEPENDENCIES(specfem sem serial both)
ENDIF(ISPARALLEL)

# ------------------------ END SPECFEM ------------------------



# Build specfemx
#TARGET_LINK_LIBRARIES(specfem ${semall})
#SET_TARGET_PROPERTIES(specfem PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

IF(ISPARALLEL)
  SET_PROPERTY(SOURCE specfemx.F90 PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE specfemx.F90 PROPERTY INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
  SET(LIBS specfem both parsem sem ${MPI_Fortran_LIBRARIES})
  SET(mainprog pspecfemx)
ELSE(ISPARALLEL)
  SET(LIBS specfem both serial sem)
  SET(mainprog specfemx)
ENDIF(ISPARALLEL)

ADD_EXECUTABLE(${mainprog}
              specfemx.F90)
TARGET_LINK_LIBRARIES(${mainprog} ${LIBS})

# This is necessary 
IF(USE_PETSC)
  SET_PROPERTY(SOURCE specfemx.F90 APPEND PROPERTY COMPILE_DEFINITIONS USE_PETSC)
  IF(ISCOMPLEX)
    SET_PROPERTY(SOURCE specfemx.F90 APPEND PROPERTY COMPILE_DEFINITIONS USE_COMPLEX)
   ENDIF(ISCOMPLEX)
ENDIF(USE_PETSC)
IF(USE_PETSC)
  TARGET_INCLUDE_DIRECTORIES(both PUBLIC "${PETSC_DIR}/include" "${PETSC_ARCH}/include")
  TARGET_LINK_LIBRARIES(both ${PETSC_ARCH}/lib/libpetsc.so)
ENDIF(USE_PETSC)




# Build partmesh program
# partmesh program partitions the mesh for the parallel processing
IF(ISPARTMESH)
  SET(srcpartmesh
     global.f90
     shared.f90
     string_library.f90
     serial_library.f90
     math_library.f90
     element.f90
     infinite_layer.f90
     partmesh_library.f90
     recreate_fault.f90
     read_input.f90
     partmesh_scotch.f90
     partmesh.f90
  )
  ADD_LIBRARY(partscotch ${srcpartmesh})
  TARGET_INCLUDE_DIRECTORIES(partscotch PUBLIC
${SCOTCH_INCLUDE_PATH})
  ADD_EXECUTABLE(partmesh
              partmesh.f90)
  TARGET_LINK_LIBRARIES(partmesh partscotch
${SCOTCH_LIBRARY_PATH}/libscotch.a ${SCOTCH_LIBRARY_PATH}/libscotcherr.a
-lpthread)
ENDIF(ISPARTMESH)
