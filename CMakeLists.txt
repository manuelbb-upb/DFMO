cmake_minimum_required(VERSION 3.12)
                                        # Added in version 3.12: Object libraries can be linked to with target_link_libraries().
project(dfmo
        DESCRIPTION "Fork of DFMO solver."
        LANGUAGES Fortran
)

enable_language(Fortran)
# Set the fortran compiler name
# set(CMAKE_Fortran_COMPILER "gfortran" CACHE INTERNAL "")

# Define flags that are used globally
set(CMAKE_Fortran_FLAGS "-O3 -g" CACHE INTERNAL "unconditional global flags")

if(WIN32)
	set(FEXT_EXE ".exe" CACHE INTERNAL "File extension for Windows executable.")
	set(FEXT_SHARED ".dll" CACHE INTERNAL "File extension for Windows shared library.")
endif()
if(UNIX)
	set(FEXT_EXE ".bin" CACHE INTERNAL "File extension for Linux executable.")
	set(FEXT_SHARED ".so" CACHE INTERNAL "File extension for Linux shared library.") 
endif()

# Define flags based on target
set(LIB_FLAGS_EXE "" CACHE INTERNAL "Flags for standalone executable target objects")
set(LIB_FLAGS_SHARED "${LIB_FLAGS_EXE}" "-fPIC" CACHE INTERNAL "Flags for shared library target objects")

set(FLAGS_EXE "" CACHE INTERNAL "Flags for standalone executable target")
set(FLAGS_SHARED "${FLAGS_EXE}" "-shared" CACHE INTERNAL "Flags for shared library target")

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
function(add_fortran_target target_name lib_flags flags problem_src)
  # Target-specific names of module libraries
  set(PROBLEM_LIB_NAME "problem_${target_name}")
  set(OBJECTS_LIB_NAME "objs_${target_name}")

  # Add common objects and set flags
  add_library(${OBJECTS_LIB_NAME} "${COMMON_OBJS}")
  target_compile_options(${OBJECTS_LIB_NAME} PRIVATE ${lib_flags})

  # Add problem module, set flags, and make it a dependency of OBJECTS_LIB_NAME
  add_library(${PROBLEM_LIB_NAME} OBJECT ${problem_src})
  target_compile_options(${PROBLEM_LIB_NAME} PRIVATE ${lib_flags})
  target_link_libraries(${OBJECTS_LIB_NAME} PUBLIC ${PROBLEM_LIB_NAME})

  # Define executable, to be built from the object files
  add_executable(${target_name} $<TARGET_OBJECTS:${OBJECTS_LIB_NAME}> $<TARGET_OBJECTS:${PROBLEM_LIB_NAME}>)
  target_link_options(${target_name} PRIVATE ${flags})
endfunction()

# Create executable target
add_fortran_target(exe "${LIB_FLAGS_EXE}" "${FLAGS_EXE}" problem_shim.f90)
# Create shared target
add_fortran_target(shared "${LIB_FLAGS_SHARED}" "${FLAGS_SHARED}" problem_abstract.f03)

# Set output names for targets
set_target_properties(exe PROPERTIES OUTPUT_NAME "multiobj")
set_target_properties(exe PROPERTIES SUFFIX "${FEXT_EXE}")
set_target_properties(shared PROPERTIES OUTPUT_NAME "multiobj")
set_target_properties(shared PROPERTIES SUFFIX "${FEXT_SHARED}")
