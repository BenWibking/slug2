
### What is SLUG ###

SLUG is the Stochastically Lighting Up Galaxies Code. SLUG is a stellar population synthesis code, and in many respects is similar to other commonly-used SPS codes such as [starburst99](http://www.stsci.edu/science/starburst99/docs/default.htm) and [FSPS](https://code.google.com/p/fsps/). Given an input star formation history, stellar initial mass function, and a set of evolutionary tracks and stellar atmospheres, all of these codes can predict the spectra and photometric properties of a stellar population. The main difference between SLUG and conventional SPS codes is that, instead of the usual approach of assuming that all stellar masses and ages are fully populated, SLUG is capable of stochastically sampling from the stellar initial mass function and age distribution, and thereby predicting not just the mean spectrum and photometry, but also the full distribution of spectra and photometry that results from stochastic sampling. This capability is critical in the regime of low star formation rates and total stellar masses, where finite sampling can lead to a distribution of properties that is extremely broad, and the mean values produced by other SPS codes are therefore of limited predictive power.

### Version history ###

This repository contains SLUG v2. SLUG v1 is available at <https://sites.google.com/site/runslug/>. The code here is a complete re-implementation of SLUG, with a number of major improvements.

* The ability to predict full spectra as well as photometry
* A vastly larger array of photometric filters, with the ability to add more without needing to alter the source code; there is also wide range of choices for photometric system
* The ability to handle a very wide range of functional forms for the initial mass function, cluster mass function, cluster lifetime function, and star formation history
* Many sampling methods available for mass-limited sampling
* The ability to simulate simple stellar populations with finite masses, as well as full star formation histories
* More accurate numerical methods for interpolation and integration
* A cleaner, easier-to-use control file structure
* Python helper routines for parsing output (replacing the older IDL routines used by SLUG v1)

SLUG v1 will no longer be maintained, so all users are encouraged to migrate to SLUG v2 as soon as possible.

### Layout of the repository ###

The SLUG repository contains several subdirectories.

* *bin*: this contains the slug source code when it is compiled
* *doc*: this contains documentation for the code
* *lib*: this contains the data files required by the code
    - *lib/atmospheres*: contains stellar atmosphere models
    - *lib/clf*: contains predefined star cluster lifetime functions
    - *lib/cmf*: contains predefined star cluster mass fuctions
    - *lib/filters*: contains predefined photometric filters
    - *lib/imf*: contains predefined initial mass functions
    - *lib/tracks*: contains stellar evolution tracks
* *output*: this is the default location for writing output
* *param*: this is the default location for parameter files that control slug simulations
* *slugPy*: this contains the python module slugPy, which includes routines for parsing simulation outputs
* *src*: this contains the main source code for SLUG
* *test*: this contains python test scripts
    - *test/aux*: contains data files for the test scripts

### Installing and configuring ###

#### Dependencies ####

* The SLUG source code requires the [GNU scientific library](http://www.gnu.org/software/gsl/)
* The SLUG source code requires the [BOOST c++ libraries](http://www.boost.org/)
* The Python parsing routines require numPy and sciPy

#### Compiling ####

SLUG comes with a Makefile in the main directory, and it shoudl be possible to build the code in most environment simply by typing "make". Compliation requires that the GSL and BOOST header files be included in the include path, and that their shared object libraries be included in the link path.

#### Running ####

Once the code is compiled, running a SLUG simulation is fairly straightforward. The required steps are:

1. Set the environment variable SLUG_DIR to the directory where you have installed SLUG. This is not strictly necessary, but it avoids the need for you to specify manually the locations of all the data files the SLUG needs.
2. Create a parameter file, or edit one of the provided examples, to set up the desired simulation. The files *param/example.param* and *param/example_cluster.param* can provide useful starting points for creating your own parameter files.
3. Run the simulation with the command line `./bin/slug param/mysimulation.param`.
4. Once the simulation is done, you can examine the output by eye (if you chose to write the output in ASCII format) or, more likely, read the data using the provided python routines. The two basic output reading routines can be invoked by doing the following in a python program of the command prompt of an interactive session:

`from slugPy import *`

`int_data = read_integrated('RUN_NAME')`

`cluster_data = read_cluster('RUN_NAME')`

### Questions, bugs, and getting involved ###

If you have questions about SLUG, have discovered any bugs, or want to contribute to ongoing development, please contact [Mark Krumholz](https://sites.google.com/a/ucsc.edu/krumholz/), mkrumhol@ucsc.edu.

### Acknowledgements ###

SLUG was the product of many people's work. The core code was originally written by Robert da Silva and Michele Fumagalli, with contributions from Jonathan Parra. The reimplementation and python parsing routines were written by Mark Krumholz. SLUG benefitted enormously from the methods and library of track and atmosphere models curated by Claus Leitherer as part of [starburst99](http://www.stsci.edu/science/starburst99/docs/default.htm). The extensive set of filters available is thanks to Charlie Conroy, who maintains the filter data as part of [FSPS](https://code.google.com/p/fsps/).