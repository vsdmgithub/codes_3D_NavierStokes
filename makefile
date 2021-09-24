# MAKEFILE FOR EULER

# COMPILER
cc=gfortran

# LIBRARY LOCATION
#lb_lc= -I/home/sugan/fftw/include
lb_lc= -I/home/coco/fftw/include

# LIBRARY FILE
lb_fftw=-L/home/coco/fftw/lib -lfftw3 -lm
#lb_fftw=-L/home/sugan/fftw/lib -lfftw3 -lm

# PROGRAM
program=NSE_code.f90

# MODULES
timer_mod            =modules-secondary/system_timer.f90
fft_mod              =modules-secondary/system_fftw.f90
vtr_mod              =modules-secondary/system_VTR.f90
vtk_mod              =modules-secondary/system_VTK.f90
constants_mod        =modules-secondary/system_constants.f90
auxilaries_mod       =modules-secondary/system_auxilaries.f90
basicvariables_mod   =modules-primary/system_basicvariables.f90
basicdeclar_mod      =modules-primary/system_basicdeclaration.f90
advvariables_mod     =modules-primary/system_advvariables.f90
advdeclaration_mod   =modules-primary/system_advdeclaration.f90
initialcondition_mod =modules-primary/system_initialcondition.f90
basicfunctions_mod   =modules-primary/system_basicfunctions.f90
advfunctions_mod     =modules-primary/system_advfunctions.f90
solver_mod          =modules-secondary/system_solver.f90
test_mod             =modules-secondary/system_test.f90
basicoutput_mod      =modules-primary/system_basicoutput.f90
advoutput_mod        =modules-primary/system_advoutput.f90
pvdoutput_mod        =modules-primary/system_pvdoutput.f90
main_mod             =modules-primary/system_main.f90

# OBJECTS
obj=system_timer.o\
	system_fftw.o\
	system_VTR.o\
	system_VTK.o\
	system_constants.o\
	system_auxilaries.o\
	system_basicvariables.o\
	system_advvariables.o\
	system_advdeclaration.o\
	system_initialcondition.o\
	system_basicfunctions.o\
	system_basicdeclaration.o\
	system_advfunctions.o\
	system_solver.o\
	system_basicoutput.o\
	system_advoutput.o\
	system_test.o\
	system_pvdoutput.o\
	system_main.o

# EXECUTABLE
run=./ex

# CLEAN COMMANDS
rmex=rm ex

mkcl=make cl
#----------------------------end-------

# MAKEFILE
# ---------------------------start-----
ex:$(ob)
	$(cc) $(lb_lc) -c $(fft_mod) $(lb_fftw)
	$(cc) -c $(timer_mod)
	$(cc) -c $(vtr_mod)
	$(cc) -c $(vtk_mod)
	$(cc) -c $(constants_mod)
	$(cc) -c $(auxilaries_mod)
	$(cc) -c $(basicvariables_mod)
	$(cc) -c $(basicdeclar_mod)
	$(cc) -c $(initialcondition_mod)
	$(cc) -c $(basicoutput_mod)
	$(cc) -c $(basicfunctions_mod)
	$(cc) -c $(solver_mod)
	$(cc) -c $(test_mod)
	$(cc) -c $(advvariables_mod)
	$(cc) -c $(advdeclaration_mod)
	$(cc) -c $(advoutput_mod)
	$(cc) -c $(pvdoutput_mod)
	$(cc) -c $(advfunctions_mod)
	$(cc) -c $(main_mod)
	$(cc) $(lb_lc) $(program) $(obj) $(lb_fftw) -o ex
	$(mkcl)
	$(run)

#----------------------------end-------

# CLEANING
# ---------------------------start-----
clean:
	rm ex
	clear
cl:
	rm *.mod
	rm *.o
#----------------------------end-------
