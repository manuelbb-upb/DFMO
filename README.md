# DFMO24
## Custom Fork of »DFMO« Multi-Objective Optimizer

### About
This fork of [DFMO](https://github.com/DerivativeFreeLibrary/DFMO) contains slight 
modifications of the original Fortran code to enable compilation to a shared library.
The shared library can be used to interactively call DFMO from other programming 
languages.
We use it in [DFMOWrapper.jl](https://github.com/manuelbb-upb/DFMOWrapper.jl) to
solve problems in the Julia language.

#### Details

Due to the nature of Fortran, an ahead-of-time compiled language, problems are hard-coded 
in the original source -- in the file `problem.f90`.
Interacting with dynamic languages is thus inhibited.
E.g., we might want to define a multi-objective problem (or many such problems) in Julia and 
call DFMO without lossy translation wizardry and recompilation.
Luckily, Fortran, like C, can be compiled to shared libraries. Shared libraries can be loaded
by other programs at runtime.
Even more luckily, Julia has built-in support for loading Fortran shared library objects!

These are the main changes compared to the original source:
* Besides compilation to a static executable, we support compilation to a shared library.
* The build system has been changed from `make` to `cmake`. Details on manual compilation  
  can be found below.
* Windows support: I have never tried to compile Fortran on Windows, so the original code
  could just easily compile on Windows, as well.  
  For our fork, I setup `cmake` to be able to cross-compile to Windows from Linux, thanks 
  to the [mingw project](https://www.mingw-w64.org/).
* The main optimization routine, formerly a Fortran `PROGRAM`, is now wrapped as a `SUBROUTINE` stored in `opt_multiobj.f90`.  
  This subroutine `opt_multiobj` is available in the shared library and takes the algorithm settings as its arguments.  
  The file `main.f90` is kept for compilation to static executables and contains a `PROGRAM` calling `opt_multiobj` with
  default settings.
* The (optimization) problem is now provided in a Fortran `MODULE`. 
  Depending on the compilation mode, different files are relevant: 
  - Static executable: The file `problem.f90` is kept for backwards compatibility, e.g. for use with TESTMO problems.
    It is, however, not directly included/required by `main.f90`.
    Rather, there is a shim module defined in `problem_shim.f90` that literally includes `problem.f90`.
    The main program is linked with the shim module.
  - Shared library: The problem module is defined in `problem_abstract.f03`.
    Functions that previously defined an optimization problem, such as `setdim` or `functs` are allocated as procedure
    pointers. There are interface definitions for all the relevant procedures.
    Furthermore, to actually enable passing callbacks from other programs (without changing function call signatures in 
    the rest of the code), there are setter functions. 
    For example, `set_setdim_ptr` is available in the shared library and takes a single argument, the pointer to a callback
    function, that is then assigned to the corresponding internal Fortran procedure pointer.
* Automatic compilation using Github Workflows. Each push to `main` triggers a build of shared libraries, and everytime a `tag` 
  is pushed (to any branch), a release is created and the build results are attached.

### Usage

* To see how problems are defined and what the algorithm setting parameters mean, please refer to `README.original.txt`.
* To use DFMO from within Julia, see [DFMOWrapper.jl](https://github.com/manuelbb-upb/DFMOWrapper.jl). 

#### (Manual) Compilation
Install `cmake` and `gfortran` to compile for Linux.
To cross-compile to Windows, install the correct mingw version of `gfortran`.
Depending on your OS, it might be called `x86_64-w64-mingw32-gfortran` or similar.

To compile the static executable, `cd` into this directory.  
Then create the output folder with `mkdir build` and change into it, `cd build`.
Configure the cmake environment: `cmake ..`, or `cmake -DCMAKE_TOOLCHAIN_FILE=../TC-mingw.cmake ..` for Windows.
Finally, compile with `cmake --build . --target "exe"`.

If I am not mistaken, it should be possible to build both targets, "shared" and "exe", within the same directory, and without reconfiguration. 
