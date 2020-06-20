#compiler options for debugging
c=gfortran -fimplicit-none -fbounds-check -fbacktrace -g -O0 -Wunused -finit-real=nan
# c=gfortran -fimplicit-none -fbounds-check -fbacktrace -g -O0 -Wunused
#compiler options for fast execution
#c=gfortran -O3

all_objs:= types.o parameters.o functions.o integ.o linalg.o simpleinteg.o

all: $(all_objs) main.o
	$c -g -o methods $(all_objs) main.o



types.o: types.f90
	$c -c types.f90

parameters.o: types.o parameters.f90
	$c -c parameters.f90

functions.o: types.o parameters.o functions.f90
	$c -c functions.f90
# here we specifying  rules for created module (integ)
# integ.o is compiled object then comes objects which inges.o dependeds on.
# has to depend on the source file next line is the recipy how to create the line.
# compiler specifications. -c - it will work with source file. and then the na,e of the source filei

integ.o: types.o functions.o linalg.o integ.f90
	$c -c integ.f90

simpleinteg.o: types.o functions.o linalg.o simpleinteg.f90
	$c -c simpleinteg.f90

linalg.o: types.o linalg.f90
	$c -c linalg.f90

main.o: $(all_objs) main.f90
	$c -c main.f90

clean:
	rm -rf *.o methods
