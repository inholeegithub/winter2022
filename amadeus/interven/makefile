cmpl=ifort
segl=ifort
#OPT = -fast
#OPT = -check all
OPT = 
OPT = -CB -check all -warn unused
OPT = -CB -check all -warn interface -assume realloc_lhs
FFLAGS = -c ${OPT}
FFLAGSc = -c 

FILES= fixedatoms.o strings.o 

/home/ihlee/csa_vasp/fixedatoms.x: $(FILES)	 
	$(segl) -o /home/ihlee/csa_vasp/fixedatoms.x $(FILES)
fixedatoms.o:fixedatoms.f90 strings.o
	$(cmpl) $(FFLAGS) fixedatoms.f90
strings.o:strings.f90
	$(cmpl) $(FFLAGS) strings.f90

clean:
	rm -f *.x *.o *.mod *.M core*  *.ps  *.ps~  *__genmod.f90  *~
touch:
	touch *.f90 *.i makefile ; chmod 600 *.f90 *.i makefile ; ls -l *.f90 *.i makefile
rmo:
	rm -f *.o *.mod *.M core*  *.ps~  *__genmod.f90  *~
lsl:
	ls -l *.f90 makefile *.i
a2ps:
	a2ps -o fixedatoms.ps fixedatoms.f90 --prologue=color 
