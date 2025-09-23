# Fourth Term Computation

## What this code is about

The C++ program **`fourth.cpp`** computes the **fourth term** in equation (3.55):

$$
\sum_{k=0}^{\infty} (-1)^k \mu^{-(2k+2)} \beta^k =
\sum_{k=0}^{\lfloor d/2\rfloor - 1} (-1)^k \beta^k \bigl(A_{2k} + B_{2k} + C_{2k}\bigr) +
\sum_{k=\lfloor d/2\rfloor}^{\infty} (-1)^k \beta^k D_{2k},
$$

where

$$
D_k =
\sum_{m=0}^{d} c_m m!
\sum_{l=0}^{m} \frac{(-1)^l}{(l!)^2 (m - l)!}
\mathrm{FP}\int_{0}^{\infty} \frac{e^{-x/2}}{ x^{k + \nu + 1 - l}} dx,
$$

and

$$
\mathrm{FP}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l} dx
= \frac{(-1)^{k-l+1}\left(\frac{1}{2}\right)^{k-l+1+\nu}
  \pi}{\Gamma(k - l + 1 + \nu)\sin(\pi\nu)}.
$$ 

The program requires $d + 1$ coefficients $c_m$ as inputs,  
read from **`Constant.txt`**, and outputs the computed values of the series for  
$\beta \in \{10^{-2}, \dots, 10^{25}, 0.2, 4\}$ to **`FOURTH.txt`**.

In the applications considered, the term (1) is typically computed up to `k = 2d`.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `fourth.cpp`      | Main C++ source computing the fourth term in equation (3.55).           |
| `Constant.txt`    | Input data: the d + 1 coefficients $c_m$.                             |
| `FOURTH.txt`      | Output file with computed series values.                                 |
| `CMakeLists.txt`  | Modern build configuration using CMake (MPI + Boost + MPFR).            |
| `compile.job`     | SLURM script to build the project on an HPC cluster.                    |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |

---

## Building and Running
To build and run the program execute shell `run.sh` script.
