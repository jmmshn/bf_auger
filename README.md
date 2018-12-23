# Brutef-force Auger Calculations (Beta)

While I have tried to make the code more user friendly, this code is still very much a Kludge.
If you have any touble getting things to work, or if you have any questions, do not hesitate to ask.

## Direct Auger recombination

- Obtain the eigenenergies and wavefunctions on a irreducible wedge on a dense *k*-point grid
- The user sets carrier concentration and the program will apply a windowing routine to obtain the Fermi-level
- Once the Fermi-level is calculated, we compute the all momentum conserving 4-body interactions (considering Umklapp) to get all of the *k*-points outside of the IBZ that need to be calculated. 
- Compute the eigenenergies and wavefunctions at the new *k*-points
- Finally compute the Auger recombination rate as a function of the band gap.

## Phonon-assisted Auger recombination

- Since the phonon-assisted process does not require momentum conservation (the phonon can provide arbitrary amount of )

## Getting Started

Clone the repository to a compute with FFTW 

```
git clone
```

### Prerequisites

Make sure that cray-FFTW is installed on the cluster you are running on.
And the module is loaded.

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

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Original code by Daniel Steiauf, 


# bf_auger
Brute force Auger calculations


# Structure of code
```
├── README.md
├── example
└── src
    ├── Makefile
    ├── base.F90
    ├── bf_auger.F90
    ├── bf_auger_sub_hhe.F90
    ├── check_kind.f90
    ├── energies.F90
    ├── inpfile.f90
    ├── kgrids.f90
    ├── klists.F90
    ├── make.depend
    ├── me_direct.F90
    ├── me_fft_direct.F90
    └── me_test.f90
```
