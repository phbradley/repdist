# repdist
Software for comparing T cell receptor (TCR) repertoires

This repository contains C++ source code for computing Repdist and TStar repertoire comparison scores for human and mouse TCR repertoires. Usage examples can be found in the shell script: `test/run.bash`

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

# EXAMPLES

Here is an example command line that compares two human beta-chain repertoires from the `test/` directory:

```
./bin/repdist --outprefix test_ --tcrs_file1 test/human_beta_tcrs_file1.txt --tcrs_file2 test/human_beta_tcrs_file2.txt --organism human --database ./db
```

And this is what the input file format looks like (simple comma-separated: V-gene,CDR3; one TCR per line; no header):

```
repdist$ head -n3 test/human_beta_tcrs_file1.txt 
TRBV10-1*01,CASSSGTANTEAFF
TRBV10-2*01,CASLDRGVHEQYF
TRBV10-2*01,CASQDRGGGTQYF
```

The `repdist` executable generates `.tsv` (tab-separated values) output files with the repdist scores and the per-TCR TStar scores, all prefixed with the string passed in via the `--outprefix` command line argument:

```
repdist$ ls test_* 
test_file1_NNdist_and_TStar_scores.tsv
test_file2_NNdist_and_TStar_scores.tsv
test_repdist_scores.tsv
```
