cmake_minimum_required(VERSION 3.12)
                                        # Added in version 3.12: Object libraries can be linked to with target_link_libraries().
project(dfmo
        DESCRIPTION "Fork of DFMO solver."
        LANGUAGES Fortran
)

enable_language(Fortran)
# Set the fortran compiler name
set(CMAKE_Fortran_COMPILER "gfortran")
# Define flags that are used globally
set(CMAKE_Fortran_FLAGS, "-O3" "-g")

# Define flags based on target
set(FLAGS_EXE "" CACHE STRING "Flags for standalone executable target")
set(FLAGS_SHARED "${FLAGS_EXE}" "-fPIC" CACHE STRING "Flags for shared library target")

# Define common source files
set(COMMON_OBJS
    main.f90
    opt_multiobj.f90
    DFMO.f90
    subroutines_DFMO.f90
    halton.f90
    sobol.f90
    qsortd.f90
    modules_DFMO.f90    # depends on `problem_shim.f90` or `problem_abstract.f03`
    CACHE INTERNAL "")

# Function to set up different targets with their respective flags and dependencies
function(add_fortran_target target_name flags problem_src)
  # Target-specific names of module libraries
  set(PROBLEM_LIB_NAME "problem_${target_name}")
  set(OBJECTS_LIB_NAME "objs_${target_name}")
  # Target-specific module flags
  set(LIB_FLAGS ${flags} "-c")

  # Add common objects and set flags
  add_library(${OBJECTS_LIB_NAME} "${COMMON_OBJS}")
  target_compile_options(${OBJECTS_LIB_NAME} PRIVATE ${LIB_FLAGS})

  # Add problem module, set flags, and make it a dependency of OBJECTS_LIB_NAME
  add_library(${PROBLEM_LIB_NAME} OBJECT ${problem_src})
  target_compile_options(${PROBLEM_LIB_NAME} PRIVATE ${LIB_FLAGS})
  target_link_libraries(${OBJECTS_LIB_NAME} PUBLIC ${PROBLEM_LIB_NAME})

  # Define executable, to be built from the object files
  add_executable(${target_name} $<TARGET_OBJECTS:${OBJECTS_LIB_NAME}> $<TARGET_OBJECTS:${PROBLEM_LIB_NAME}>)
  target_compile_options(${target_name} PRIVATE ${flags})
endfunction()

# Create executable target
add_fortran_target(exe "${FLAGS_EXE}" problem_shim.f90)
# Create shared target
add_fortran_target(shared "${FLAGS_SHARED}" problem_abstract.f03)

# Set output names for targets
set_target_properties(exe PROPERTIES OUTPUT_NAME multiobj.bin)
set_target_properties(shared PROPERTIES OUTPUT_NAME multiobj.so)
