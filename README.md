# repdist
Software for comparing T cell receptor (TCR) repertoires

This repository contains C++ source code for computing Repdist and TStar repertoire comparison scores.
Usage examples can be found in the shell script: `tests/run.bash`

# REQUIREMENTS

This software depends on header files included with the BOOST C++ library.
You can download the library [here](https://www.boost.org/users/download/).

# COMPILING

Edit the "BOOSTDIR" line in the Makefile (in the top directory) to point to the location where your BOOST download is installed. Then type `make` (while in the top, repdist/ directory). The binary executable files will be placed in the `bin/` directory.

# THANKS

We are using the [TCLAP](http://tclap.sourceforge.net/) header library for parsing command line arguments. As suggested by the TCLAP docs, we have included the header files within this repository for ease of compiling. Please see the author and license information in `include/tclap/`.

# TESTING

There is a simple bash script that tests the `repdist` executable. To run it:

```
cd test/
./run.bash
```

