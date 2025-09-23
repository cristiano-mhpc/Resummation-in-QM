$$\newcommand{\bbint}[2]{\;\backslash\!\!\!\!\backslash\!\!\!\!\!\int_{#1}^{#2}}$$ 
# Second Term Computation (MPI + Boost.MPFR)

## What this code is about

The C++ code **`second.cpp`** computes the **second term** in equation (3.10):

$$
\sum_{k=0}^{\infty} (-1)^k \mu^{-(k+1)} \beta^k=
\sum_{k=0}^{d} (-1)^k \beta^k (A_k + B_k + C_k)
+
\sum_{k=d+1}^{\infty} (-1)^k \beta^k D_k.
\tag{1}
$$

where

$$
B_k =
\sum_{m=k+1}^{d} c_m m!
\sum_{l=0}^{k}
\frac{(-1)^l}{(l!)^2 (m - l)!}
\bbint{0}{\infty} \frac{e^{-x/2}}{ x^{k + \nu + 1 - l}}\, dx,
\tag{2}
$$

and the integral is

$$
\bbint{0}{\infty} e^{-x/2} x^{k+\nu+1-l}\, dx=
\frac{(-1)^{k-l+1}
\left(\frac{1}{2}\right)^{k-l+1+\nu}
\pi}{\,\Gamma(k - l + 1 + \nu)\,\sin(\pi\nu)}.
\tag{3}
$$

The code requires the $d + 1$ numbers $c_m$ as inputs,  read from the file **`Constants.txt`**.  It outputs the computed values for various $\beta$ to **`SECOND.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `second.cpp`      | C++ source computing the second term in equation (3.10).                |
| `Constants.txt`   | Input data: the $d + 1$ coefficients $c_m$.                              |
| `SECOND.txt`      | Output file with computed series values.                                 |
| `CMakeLists.txt`  | Modern build configuration (MPI + Boost + MPFR + GMP).                  |
| `compile.job`     | SLURM script to compile the code on an HPC cluster.                     |
| `together.job`    | SLURM script to run the executable on an HPC cluster.                   |
| `mpfr.sh`         | Legacy local shell script to compile and run on Ubuntu 22.04.           |

---

## Building and Running

To build and run the code, execute the script `run.sh`. 
