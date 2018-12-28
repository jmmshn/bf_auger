# The compiler
FC = ftn
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -fast -DMPI
INCLUDES = -I/opt/cray/pe/fftw/3.3.6.2/haswell/include
LFLAGS = -L/opt/cray/pe/fftw/3.3.6.2/haswell/lib
LIBS = -lfftw3

all: eeh hhe

# library flags go at the end of main progam
eeh : AUGCHG = -DEEH
eeh : tableio inpfile base me_fft calculate_auger_phonon
	$(FC) $(INCLUDES) -o calc_auger_EEH.mpi tableio.o inpfile.o base.o me_fft.o calculate_auger_phonon.o $(LFLAGS) $(LIBS)

hhe : AUGCHG = -DHHE
hhe : tableio inpfile base me_fft calculate_auger_phonon
	$(FC) $(INCLUDES) -o calc_auger_HHE.mpi tableio.o inpfile.o base.o me_fft.o calculate_auger_phonon.o $(LFLAGS) $(LIBS)

me_fft : me_fft.F90
	$(FC) ${FCFLAGS} ${AUGCHG} -c $<
tableio : tableio.f90
	$(FC) ${FCFLAGS} ${AUGCHG} -c $<
calculate_auger_phonon : calculate_auger_phonon.F90
	$(FC) ${FCFLAGS} ${AUGCHG} -c $<
inpfile : inpfile.f90
	$(FC) ${FCFLAGS} ${AUGCHG} -c $<
base    : base.F90 
	$(FC) ${FCFLAGS} ${AUGCHG} -c $<

clean:
	rm *.o *.mod
