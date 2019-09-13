# repdist
Software for comparing T cell receptor (TCR) repertoires

This repository contains C++ source code for computing Repdist and TStar repertoire comparison scores.
Usage examples can be found in the shell script: `test/run.bash`

# COMPILING

From the command line type `make` while in the top, `repdist` directory. The binary executable files will be placed in the `bin/` directory. You may need to edit the `Makefile` to point to your C++ compiler.

Alternatively, you can just run this command from the top directory:

`g++ -O3 -std=c++11 -Wall -I ./include/  -o bin/repdist src/repdist.cc`

changing `g++` to your C++ compiler if necessary.

# THANKS

We are using the [TCLAP](http://tclap.sourceforge.net/) header library for parsing command line arguments. As suggested by the TCLAP docs, we have included the header files within this repository for ease of compiling. Please see the author and license information in `include/tclap/`.

# TESTING

There is a simple bash script that tests the `repdist` executable. To run it:

```
cd test/
./run.bash
```

