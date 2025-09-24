# Fourth Term Computation (MPI + Boost.MPFR)

## What this code is about

The C++ code **`fourth.cpp`** computes the **fourth term** in equation (3.10):

$$\sum_{k=0}^{\infty} (-1)^k \mu^{-(k+1)} \beta^k=\sum_{k=0}^{d} (-1)^k \beta^k (A_k + B_k + C_k)+\sum_{k=d+1}^{\infty} (-1)^k \beta^k D_k,$$

where

$$D_k =\sum_{m=0}^{d} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m - l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l} dx,$$

and

$$\mathrm{FP}\int_{0}^{\infty} \frac{e^{-x/2}}{ x^{k+\nu+1-l}} dx=\frac{(-1)^{k-l+1}\left(\frac{1}{2}\right)^{k-l+1+\nu}\pi}{ \Gamma(k - l + 1 + \nu)\sin(\pi \nu)}.$$

The program reads the $d + 1$ coefficients $c_m$ from **`Constants.txt`**,and outputs the computed values of the fourth-term series to **`FOURTH.txt`**  for a range of $\beta$ values.In practice, the term (1) is typically computed up to $k = 2d$.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `fourth.cpp`      | C++ source computing the fourth term $D_k$ in equation (3.10).         |
| `Constants.txt`   | Input file with the $d + 1$ coefficients $c_m$.                         |
| `FOURTH.txt`      | Output file with computed series values.                                |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM script to compile the program on an HPC cluster.                  |
| `together.job`    | SLURM script to run the executable on an HPC cluster.                   |
| `run.sh`         | For running and building the program.                 |

---

## Building and Running

### 1. Local Ubuntu build with CMake
To build and run the program, execute the `run.sh` 
