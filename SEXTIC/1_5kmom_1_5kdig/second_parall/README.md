$$\newcommand{\bbint}[2]{\;\backslash\!\!\!\!\backslash\!\!\!\!\!\int_{#1}^{#2}}$$ 

# Second Term Computation

## What this code is about

The C++ program **`second.cpp`** computes the **second term** in equation (3.55):

$$
\sum_{k=0}^{\infty} (-1)^k \mu^{-(2k+2)} \beta^k =
\sum_{k=0}^{\lfloor d/2\rfloor - 1} (-1)^k \beta^k \bigl(A_{2k} + B_{2k} + C_{2k}\bigr) +
\sum_{k=\lfloor d/2\rfloor}^{\infty} (-1)^k \beta^k D_{2k},
$$

where $\lfloor x \rfloor$ is the floor function and  
the $B_k$ term is given by

$$
B_k =
\sum_{m=k+1}^{d} c_m m! 
\sum_{l=0}^{k} \frac{(-1)^l}{(l!)^2 (m - l)!}
\int_{0}^{\infty} e^{-x/2} x^{k + \nu + 1 - l} dx,
$$

and the integral evaluates to

$$
\bbint{0}{\infty} \frac{e^{-x/2}}{ x^{k+\nu+1-l}} dx
= \frac{(-1)^{k-l+1}\left(\frac{1}{2}\right)^{k-l+1+\nu}
  \pi}{\Gamma(k - l + 1 + \nu)\,\sin(\pi\nu)}.
$$

The program requires the **d + 1 coefficients \(c_m\)** as inputs,  
read from **`Constant.txt`**, and outputs the computed series values for  
$\beta \in \{10^{-5},\dots,10^{23},0.2,4\}$ to **`SECOND.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `second.cpp`      | Main C++ source computing the second term in equation (3.55).          |
| `Constant.txt`    | Input data containing the d + 1 coefficients $c_m$.                   |
| `SECOND.txt`      | Output file with computed series values.                                |
| `CMakeLists.txt`  | Modern build configuration using CMake (MPI + Boost + MPFR).            |
| `run.sh`          | Portable build/run script: configures, builds, and launches locally or on SLURM. |
| `compile.job`     | SLURM script to build the project on an HPC cluster.                    |
| `together.job`    | SLURM script to run the executable on an HPC cluster.                   |

---

## Building and Running
To build and run the program run the shell script `run.sh`. 
