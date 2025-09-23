# Matrix Computation for P(n,m) (MPI + Boost.MPFR)

## What this code is about

The C++ code **`LU.cpp`** computes the matrix in equation (3.8):

$$
P(n,m)
= m!\,2^{n-\nu+1}
  \sum_{k=0}^{m}
    \frac{(-2)^k \Gamma(n + k - \nu + 1)}
         {(k!)^2 (m - k)!}.
\tag{1}
$$

The code writes the elements $P(n,m)$ to the file **`matrix_p.txt`**  
in a row-major fashion.

---

## Repository Layout

| File / Folder   | Purpose                                                                 |
|-----------------|-------------------------------------------------------------------------|
| `LU.cpp`        | C++ source code computing the matrix $P(n,m)$ of equation (3.8).        |
| `matrix_p.txt`  | Output file containing the matrix elements in row-major order.          |
| `CMakeLists.txt`| Modern CMake build configuration to build `LU.cpp`.                     |
| `compile.job`   | SLURM script to compile the code on an HPC cluster.                     |
| `together.job`  | SLURM script to run the compiled executable on an HPC cluster.          |
| `run.sh`       | Legacy shell script to compile and run on a local Ubuntu 22.04 machine. |

---

## Building and Running

### 1. Local Ubuntu build with CMake
To build and run the program, execute the script `run.sh`. 
