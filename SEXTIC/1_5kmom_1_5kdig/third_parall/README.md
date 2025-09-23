$$\newcommand{\bbint}[2]{\;\backslash\!\!\!\!\backslash\!\!\!\!\!\int_{#1}^{#2}}$$ 

# Third Term Computation (MPI + Boost.MPFR)

## What this code is about

The C++ program **`third.cpp`** computes the **third term** in equation (3.55):

$$
\sum_{k=0}^{\infty} (-1)^k \mu^{-(2k+2)} \beta^k =
\sum_{k=0}^{\lfloor d/2\rfloor - 1} (-1)^k \beta^k \bigl(A_{2k} + B_{2k} + C_{2k}\bigr) +
\sum_{k=\lfloor d/2\rfloor}^{\infty} (-1)^k \beta^k D_{2k},
$$ 

where $\lfloor x \rfloor$ is the floor function and the $C_k$ term is

$$
C_k =
\sum_{m=k+1}^{d} c_m m!
\sum_{l=k+1}^{m}
\frac{(-1)^l \Gamma(l - k - \nu) 2^{\,l-k-\nu}}{(l!)^2 (m - l)!}.
$$ 

The program reads the required $d + 1$ constants $c_m$ from  
**`Constant.txt`**, and outputs the computed series values for  
$\beta \in \{10^{-5},\dots,10^{23},0.2,4\}$ to **`THIRD.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `third.cpp`       | Main C++ source computing the third term in equation (3.55).            |
| `Constant.txt`    | Input data: the d + 1 coefficients $c_m$.                             |
| `THIRD.txt`       | Output file with computed series values.                                 |
| `CMakeLists.txt`  | Modern build configuration using CMake (MPI + Boost + MPFR).            |
| `compile.job`     | SLURM script to build the project on an HPC cluster.                    |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |

---

## Building and Running
To buils and run the program run the script `run.sh`. 
