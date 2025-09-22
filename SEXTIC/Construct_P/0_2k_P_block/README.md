# LU Matrix Computation (MPI + Boost.MPFR)

## What this code is about

The C++ code **`LU.cpp`** computes the matrix elements in equation (3.54) in the paper:

$$
P(n,m) = \, 2^{2n - \nu + 1}m!
\sum_{k=0}^{m}
\frac{(-2)^k \, \Gamma(2n + k - \nu + 1)}{(k!)^2 (m - k)!}.
$$

The code writes the elements $P(n,m)$ to the file `matrix_p.txt`  in a **row-major** fashion.

---


## Supporting Files

| File            | Purpose                                                                 |
|-----------------|-------------------------------------------------------------------------|
| `compile.job`   | SLURM script to compile the code on an HPC system and generate an executable. |
| `together.job`  | SLURM script to run the executable on an HPC system.                   |
| `mpfr.sh`       | Shell script to compile and run the code on an Ubuntu 22.04 local machine. |
| `LU.cpp`        | C++ source code that computes the matrix elements using MPI and Boost.MPFR. |
| `run.sh`       | Shell script to alternatively build and run the executable on using CMake.                    |
| `CMakeLists.txt` | CMake configuration file to build the project.                          |``
---

## Quick Start with Existing Script (Local)

On a local Ubuntu machine (with MPI, Boost, and MPFR installed),  
you can compile and run as follows:

```bash
chmod +x mpfr.sh
./mpfr.sh

