#compiler options for debugging
c=gfortran -fimplicit-none -fbounds-check -fbacktrace -g -O0 -Wunused -finit-real=nan
# c=gfortran -fimplicit-none -fbounds-check -fbacktrace -g -O0 -Wunused
#compiler options for fast execution
#c=gfortran -O3

all_objs:= types.o parameters.o functions.o integ.o

all: $(all_objs) main.o
	$c -g -o methods $(all_objs) main.o

types.o: types.f90
	$c -c types.f90

parameters.o: types.o parameters.f90
	$c -c parameters.f90

functions.o: types.o parameters.o functions.f90
	$c -c functions.f90

integ.o: types.o functions.o integ.f90
	$c -c integ.f90

main.o: $(all_objs) main.f90
	$c -c main.f90

clean:
	rm -rf *.o methods
