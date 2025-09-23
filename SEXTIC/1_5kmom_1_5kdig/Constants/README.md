# Computation of Expansion Coefficients $c_m$

## What this code is about


The C++ file **`LU.cpp`** solves the system of $d+1$ linear equations  
for the first $d+1$ expansion coefficients $c_m$:
$\mu_n =\sum_{m=0}^{d} c_m P(n,m)\tag{1}$
where the matrix $P(n,m)$ is defined as

$$P(n,m) = m!2^{n-\nu+1}\sum_{k=0}^{m}\frac{(-2)^k\Gamma(n + k - \nu + 1)}{(k!)^2 (m - k)!}$$
In the case of the **quartic anharmonic oscillator**, the coefficients are related by
$b_{k+1} = (-1)^k \mu_k
       = (-1)^k \int_{0}^{\infty} x^k \rho(x)\,dx,
\tag{3}$
for $k = 0, 1, \dots$, where $b_{k+1}$ are the coefficients of the  
weak-coupling perturbation expansion for the ground-state energy
$E^{(2)}(\beta) = 1 + \sum_{k=1}^{\infty} b_k \beta^k.
\tag{4}$
The program outputs the coefficients $c_m$ and related quantities needed  
in the weak-coupling expansion.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `LU.cpp`          | C++ source code that solves the linear system and computes $c_m$.     |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM script to compile the program on an HPC cluster.                  |
| `together.job`    | SLURM script to run the program on an HPC cluster.                      |

---

## Building and Running
To build and run the program, execute the shell script `run.sh`. 
