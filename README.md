# Brute-Force Auger Calculations (Beta)

While I have tried to make the code more user friendly, this code is still very much a Kludge.
If you have any touble getting things to work, or if you have any questions, do not hesitate to message me on github.

## Direct Auger recombination

- Obtain the eigenenergies and wavefunctions on a irreducible wedge on a dense *k*-point grid
- The user sets carrier concentration and the program will apply a windowing routine to obtain the Fermi-level
- Once the Fermi-level is calculated, we compute the all momentum conserving 4-body interactions (considering Umklapp) to get all of the *k*-points outside of the IBZ that need to be calculated. 
- Compute the eigenenergies and wavefunctions at the new *k*-points
- Finally compute the Auger recombination rate as a function of the band gap.

## Phonon-assisted Auger recombination

- Since the phonon-assisted process does not require momentum conservation (the phonon can provide arbitrary amount of momentum), the combinatorice will be prohibatively large if we are considering *k*-dependent filling of the of the band edges.
- We have to overcome this by assuming that all of the states are at a single point at the band-edge (or a few if they are degenerate).
- Then we compute the electron-phonon matrix elements and Coulomb matrix element
- Note that use use a patched QE to print the electron-phonon matrix elements within QE because it's important for the phase of all wavefunctions to be consistent.

### Prerequisites

Make sure that cray-FFTW is installed on the cluster you are running on and make sure the module is loaded.

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Acknowledgments

* Original code by Emmanouil Kioupakis, Daniel Steiauf, Patrick Rinke and Kris T Delaney
* 

# Structure of code
```
├── README.md
└── src
    └── direct
        ├── base.F90
        ├── bf_auger.F90
        ├── bf_auger_sub_hhe.F90
        ├── check_kind.f90
        ├── energies.F90
        ├── inpfile.f90
        ├── kgrids.f90
        ├── klists.F90
        ├── make.depend
        ├── Makefile
        ├── me_direct.F90
        ├── me_fft_direct.F90
        └── me_test.f90
```
