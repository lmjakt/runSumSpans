# Functions related to finding spans or regions from numeric data

runSumSpans using an accumulating sum to identify regions in a similar
manner to cpg_distribute written by Gos Micklem + Tim Cutts (2004) but
taking a sequence of scores, their associated positions, and a
separation penalty rather than a sequence of nucleotide characters.

## Usage

Compile the function using:

```sh
R CMD SHLIB runSumSpans.cpp
```

then from an R session:

```R
dyn.load("<path to so file>/runSumSpans.cpp")
## where path to so file is wherever you put
## runSumSpans.so after compiling.
```

Then follow examples in example.R. Note that 'example.rds' contains a
vector of numeric values obtained using code designed to quickly identify
inversions in long sequences.


**WARNING**

runSumSpans is inherently unsafe. I have not been able to get it to check the
number of arguments supplied from .Call(). Giving the wrong number of arguments
results in a segmentation fault.

Overcoming this probably requires Rcpp; or writing a package, or something.
At the moment, it seems as if R_init_runSumSpans() never gets called at all.
My understanding is that it is supposed to be called on dyn.load().

