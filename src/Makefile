# MAKEFILE for bf_auger code
# tested on CORI
FC = ftn
FCFLAGS = -DMPI # -g -O0
LIB = -L/opt/cray/pe/fftw/3.3.4.10/haswell/lib
INC = -I/opt/cray/pe/fftw/3.3.4.10/haswell/include

all : klists.x auger_eeh.x auger_hhe.x

klists.x : klists.o 
	$(FC) -o klists.x klists.o base.o energies.o kgrids.o inpfile.o \
		$(INC) $(LIB)

auger_eeh.x : bf_auger_eeh.o
	$(FC) -o auger_eeh.x me_fft_direct.o bf_auger_eeh.o base.o \
		kgrids.o inpfile.o $(INC) $(LIB)
auger_hhe.x : bf_auger_hhe.o
	$(FC) -o auger_hhe.x me_fft_direct.o bf_auger_hhe.o base.o \
		kgrids.o inpfile.o $(INC) $(LIB)
	#check_kind : kind.o 
	#me_test : kind.o inpfile.o base.o

%: %.o
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) $(INC) -c $< $(LIB) 

%.o: %.F90               
	$(FC) $(FCFLAGS) $(INC) -c $< $(LIB) 

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

include make.depend

# special rules
bf_auger_eeh.o: bf_auger.F90
	$(FC) $(FCFLAGS) -DEEH $(INC) -c -o bf_auger_eeh.o $< $(LIB) 

bf_auger_hhe.o: bf_auger.F90
	$(FC) $(FCFLAGS) $(INC) -c -o bf_auger_hhe.o $< $(LIB) 
