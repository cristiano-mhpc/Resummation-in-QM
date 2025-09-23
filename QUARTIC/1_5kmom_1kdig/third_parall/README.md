# Third Term Computation (MPI + Boost.MPFR)

## What this code is about

The C++ code **`third.cpp`** computes the **third term** in equation (3.10):

$$
\sum_{k=0}^{\infty} (-1)^k \mu^{-(k+1)} \beta^k=
\sum_{k=0}^{d} (-1)^k \beta^k (A_k + B_k + C_k)
+
\sum_{k=d+1}^{\infty} (-1)^k \beta^k D_k,
\tag{1}
$$

where

$$
C_k =
\sum_{m=k+1}^{d} c_m m!
\sum_{l=k+1}^{m}
\frac{(-1)^l \Gamma(l - k - \nu) 2^{\,l - k - \nu}}
     {(l!)^2 (m - l)!}.
\tag{2}
$$

The code requires the $d + 1$ numbers $c_m$ as inputs, read from  
**`Constants.txt`**, and writes the computed values for various $\beta$ to  
**`THIRD.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `third.cpp`       | C++ source computing the third term $C_k$ in equation (3.10).           |
| `Constants.txt`   | Input file with the $d + 1$ coefficients $c_m$.                          |
| `THIRD.txt`       | Output file containing the computed series values.                      |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM script to compile the program on an HPC cluster.                  |
| `together.job`    | SLURM script to run the program on an HPC cluster.                      |
| `mpfr.sh`         | Legacy shell script to compile and run on a local Ubuntu 22.04 machine. |

---

## Building and Running

### 1. Local Ubuntu build with CMake

Build and run by executing the script `run.sh`. 
