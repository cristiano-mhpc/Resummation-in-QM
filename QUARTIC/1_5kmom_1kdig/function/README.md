# Function Term Computation (MPI + Boost.MPFR)

## What this code is about

The C++ code **`function.cpp`** computes the **second term on the right-hand side** of equation (2.8):

$$\beta\int_{0}^{\infty}\frac{x^{-\nu} g(x)}{1 + \beta x} dx =\sum_{k=0}^{\infty} (-1)^k \mu^{-(k+1)} \beta^k+\frac{\pi g(-1/\beta) \beta^{\nu}}{\sin(\pi \nu)},$$

where

$$g(x) =e^{-x/2}\sum_{m=0}^{d} c_m m!\sum_{k=0}^{m}\frac{(-x)^k}{(k!)^2 (m - k)!}.$$

The code requires the $d + 1$ numbers $c_m$ as inputs.  
These are read from **`Constants.txt`** and the code outputs the computed values for various $\beta$ to **`FIFTH.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `function.cpp`    | C++ source computing the second term in equation (2.8).                 |
| `Constants.txt`   | Input data containing the $d + 1$ coefficients $c_m$.                   |
| `FIFTH.txt`       | Output file containing the computed results.                             |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM script to compile the code on an HPC cluster.                     |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |
| `run.sh`         | For building and running the program. |

---

## Building and Running
To build and run the code, execute the `run.sh` script. 
