# HMM for GC content emissions

To reproduce the workflow, the `Part1.R` and `Part2.R` should be run from the command line. `Part1.R` takes two arguments: the first a parameter file and the second a file of emissions.

```
Rscript Part1.R path/to/param_file.txt path/to/emission_file.txt
```

The `Part2.R` script takes a single argument - a parameter file.

```
Rscript Part2.R path/to/param_file.txt 
```

The parameter file should be formatted as follows:
- Line 1: Hidden state space S
- Line 2: Emission space V
- Line 3: Initial state distribution <img src="https://render.githubusercontent.com/render/math?math=\mu^0">
- Line 4: Transition matrix A
- Line 5: Emission matrix B

The matrices on lines 4 and 5 are filled in by row: so an entry of "0.8 0.2 0.1 0.9" corresponds to a 2x2 matrix with first row (0.8, 0.2) and second row (0.1, 0.9).
