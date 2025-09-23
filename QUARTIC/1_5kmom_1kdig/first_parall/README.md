# First Term Computation (MPI + Boost.MPFR)

## What this code is about

$$\newcommand{\bbint}[2]{\;\backslash\!\!\!\!\backslash\!\!\!\!\!\int_{#1}^{#2}}$$ 

The C++ code **`first.cpp`** computes the **first term** in equation (3.10):

$$\sum_{k=0}^{\infty}(-1)^k \mu^{-(k+1)} \beta^k=\sum_{k=0}^{d}(-1)^k \beta^k (A_k + B_k + C_k)+\sum_{k=d+1}^{\infty}(-1)^k \beta^k D_k$$

where

$$A_k =\sum_{m=0}^{k}c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m - l)!}\bbint{0}{\infty} \frac{e^{-x/2}} {x^{k+\nu+1-l}}\, dx,$$

and

$$\bbint{0}{\infty} \frac{e^{-x/2}} {x^{k+\nu+1-l} }\, dx=\frac{(-1)^{k-l+1}\left(\frac{1}{2}\right)^{k-l+1+\nu}\pi}{\Gamma(k - l + 1 + \nu)\sin(\pi \nu)}.$$

The code reads the required $d + 1$ numbers $c_m$ from **`Constants.txt`** and writes the computed values of the first-term series to **`FIRST.txt`** for a range of $\beta$ values.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `first.cpp`       | Main C++ source computing the first term in equation (3.10).            |
| `Constants.txt`   | Input data: the $d + 1$ coefficients $c_m$.                              |
| `FIRST.txt`       | Output file containing the computed first-term series values.           |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM script to compile the program on an HPC cluster.                  |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |
| `mpfr.sh`         | Legacy shell script to compile and run locally on Ubuntu 22.04.         |

---

## Building and Running

### 1. Local Ubuntu build with CMake

To build and run the program, run the shell script `run.sh`. 
