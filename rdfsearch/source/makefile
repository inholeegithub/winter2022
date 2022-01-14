cmpl=mpiifort
segl=mpiifort
#OPT = -fast
#OPT = -check all
OPT = -CB -check all -warn unused
OPT = -CB -check all -warn interface -assume realloc_lhs
OPT = -fast
OPT = -O2
FFLAGS = -c ${OPT}
FFLAGSc = -c 

FILES= invsearch.o pearson.o rcrystal.o rwycpos.o spgrp_sites.o spglattice.o strings.o rmarin.o sortnr.o numeral.o

/home/ihlee/winter2022/rdfsearch/invsearch.x: $(FILES)	 
	$(segl) -o /home/ihlee/winter2022/rdfsearch/invsearch.x $(FILES)
invsearch.o:invsearch.f90
	$(cmpl) $(FFLAGS) invsearch.f90
pearson.o:pearson.f90
	$(cmpl) $(FFLAGS) pearson.f90
rcrystal.o:rcrystal.f90  spgrp_sites.o strings.o
	$(cmpl) $(FFLAGS) rcrystal.f90
spgrp_sites.o:spgrp_sites.f90
	$(cmpl) $(FFLAGS) spgrp_sites.f90
rwycpos.o:rwycpos.f90
	$(cmpl) $(FFLAGS) rwycpos.f90
spglattice.o:spglattice.f90
	$(cmpl) $(FFLAGS) spglattice.f90
rmarin.o:rmarin.f90
	$(cmpl) $(FFLAGS) rmarin.f90
sortnr.o:sortnr.f90
	$(cmpl) $(FFLAGS) sortnr.f90
strings.o:strings.f90
	$(cmpl) $(FFLAGS) strings.f90
numeral.o:numeral.f
	$(cmpl) $(FFLAGS) numeral.f

clean:
	rm -f *.x *.o *.mod *.M core*  *.ps  *.ps~  *__genmod.f90  *~
touch:
	touch *.f90 *.i makefile ; chmod 600 *.f90 *.i makefile ; ls -l *.f90 *.i makefile
rmo:
	rm -f *.o *.mod *.M core*  *.ps~  *__genmod.f90  *~
wc:
	wc invsearch.f90 pearson.f90 rcrystal.f90 rwycpos.f90 spgrp_sites.f90 spglattice.f90
lsl:
	ls -l *.f90 makefile *.i
a2ps:
	a2ps -o invsearch.ps invsearch.f90 --prologue=color
a2psblack:
	a2ps -o invsearch.ps invsearch.f90 --medium=A4 --line-numbers=1 
ps2pdf:
	ps2pdf -sPAPERSIZE=a4 invsearch.ps 
