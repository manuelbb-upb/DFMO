x86_64-w64-mingw32-gfortran -g -shared -fPIC -o3 problem_abstract.f03 -o problem.dll

for fn in "modules_DFMO" "qsortd" "sobol" "halton" "subroutines_DFMO" "DFMO" "opt_multiobj" "main"
do 
	echo "Compiling ${fn}."
	x86_64-w64-mingw32-gfortran -c -g -O3 -shared -fPIC "${fn}.f90" -o "${fn}.dll" 
done
x86_64-w64-mingw32-gfortran -g -O3 -shared -fPIC main.dll opt_multiobj.dll DFMO.dll subroutines_DFMO.dll halton.dll sobol.dll qsortd.dll modules_DFMO.dll problem.dll -o multiobj.dll
