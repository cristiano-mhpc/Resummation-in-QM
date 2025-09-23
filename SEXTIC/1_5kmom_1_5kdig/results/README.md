# Ground-State Energy Computation (MPI + Boost.MPFR)

## What this code is about

The C++ program **`result.cpp`** computes the ground-state energy $E^{(3)}_0(\beta)$ of the **sextic anharmonic oscillator**:

$$E^{(3)}(\beta)= 1 + \sum_{k=0}^{\infty} (-1)^k \beta^k \mu^{-(2k+2)} + \Delta(\beta)$$

---

### Second term of equation (1)

$$\sum_{k=0}^{\infty} (-1)^k \mu^{-(2k+2)} \beta^k=\sum_{k=0}^{\lfloor d/2 \rfloor -1} (-1)^k \beta^k \bigl(A_{2k}+B_{2k}+C_{2k}\bigr)+\sum_{k=\lfloor d/2 \rfloor}^{\infty} (-1)^k \beta^k D_{2k}.$$

where

$$A_k =\sum_{m=0}^{k} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m-l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx,$$

$$B_k =\sum_{m=k+1}^{d} c_m m!\sum_{l=0}^{k}\frac{(-1)^l}{(l!)^2 (m-l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx,$$

$$C_k =\sum_{m=k+1}^{d} c_m m!\sum_{l=k+1}^{m}\frac{(-1)^l \Gamma(l-k-\nu) 2^{l-k-\nu}}{(l!)^2 (m-l)!},$$

$$D_k =\sum_{m=0}^{d} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m-l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx.$$

and the finite–part integral is

$$\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx= (-1)^{k-l+1}\left(\frac12\right)^{k-l+1+\nu}\pi\Gamma(k-l+1+\nu)\sin(\pi\nu).$$

Here $\nu = -\tfrac{1}{2}$.

---

### Third term of equation (1)

$$\Delta(\beta)= \frac{\pi \beta^{(1+\nu)/2}}{\sin(\pi\nu)}\Bigl[\cos\!\left(\frac{\pi\nu}{2}\right)\operatorname{Im} g\!\left(\frac{i}{\sqrt{\beta}}\right)+ \sin\!\left(\frac{\pi\nu}{2}\right)\operatorname{Re} g\!\left(\frac{i}{\sqrt{\beta}}\right)\Bigr],$$

where

$$
g(x)
= e^{-x/2}
  \sum_{m=0}^{\infty} c_m m!
  \sum_{k=0}^{m}
    \frac{(-x)^k}{(k!)^2 (m-k)!}.
\tag{9}
$$

---

## Inputs and Outputs

* Reads the precomputed series terms from  
  **`FIRST.txt`**, **`SECOND.txt`**, **`THIRD.txt`**, **`FOURTH.txt`**, and **`FIFTH.txt`**.
* Computes $E^{(3)}_0(\beta)$ for  
  $\beta \in \{10^{-2},\dots,10^{25},0.2,4\}$.
* Writes the final results to **`disectres.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `result.cpp`      | Main C++ source computing $E^{(3)}_0(\beta)$.                             |
| `FIRST.txt` … `FIFTH.txt` | Input files for the individual series terms.                     |
| `disectres.txt`   | Output file with computed ground-state energies.                         |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM script to build the project on an HPC cluster.                    |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |

---

## Building and Running
 To build and run the program, execute the script `run.sh`. 
