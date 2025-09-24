# Ground-State Energy Computation (MPI + Boost.MPFR)

## What this code is about

The C++ program computes the ground-state energy $E^{(3)}_0(\beta)$ of the **sextic anharmonic oscillator**:

$$E^{(3)}(\beta)= 1 + \sum_{k=0}^{\infty} (-1)^k \beta^k \mu^{-(2k+2)} + \Delta(\beta)$$

---

Second term in the equation above is given by 

$$\sum_{k=0}^{\infty} (-1)^k \mu^{-(2k+2)} \beta^k=\sum_{k=0}^{\lfloor d/2 \rfloor -1} (-1)^k \beta^k \bigl(A_{2k}+B_{2k}+C_{2k}\bigr)+\sum_{k=\lfloor d/2 \rfloor}^{\infty} (-1)^k \beta^k D_{2k}.$$

where

$$A_k =\sum_{m=0}^{k} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m-l)!}\mathrm{FP}\int_{0}^{\infty} \frac{e^{-x/2}}{ x^{k+\nu+1-l}}dx,$$

$$B_k =\sum_{m=k+1}^{d} c_m m!\sum_{l=0}^{k}\frac{(-1)^l}{(l!)^2 (m-l)!}\mathrm{FP}\int_{0}^{\infty} \frac{e^{-x/2} }{x^{k+\nu+1-l}}dx,$$

$$C_k =\sum_{m=k+1}^{d} c_m m!\sum_{l=k+1}^{m}\frac{(-1)^l \Gamma(l-k-\nu) 2^{l-k-\nu}}{(l!)^2 (m-l)!},$$

$$D_k =\sum_{m=0}^{d} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m-l)!}\mathrm{FP}\int_{0}^{\infty} \frac{e^{-x/2}}{ x^{k+\nu+1-l}}dx.$$

and the finite–part integral is

$$\mathrm{FP}\int_{0}^{\infty} \frac{e^{-x/2}}{ x^{k+\nu+1-l}}dx= \frac{(-1)^{k-l+1}\left(\frac12\right)^{k-l+1+\nu}\pi}{\Gamma(k-l+1+\nu)\sin(\pi\nu)}.$$

Here $\nu = -\tfrac{1}{2}$.

---

The third term is given by 

$$\Delta(\beta)= \frac{\pi \beta^{(1+\nu)/2}}{\sin(\pi\nu)}\Bigl[\cos\left(\frac{\pi\nu}{2}\right)\mathrm{Im} g\left(\frac{i}{\sqrt{\beta}}\right)+ \sin\left(\frac{\pi\nu}{2}\right)\mathrm{Re} g\left(\frac{i}{\sqrt{\beta}}\right)\Bigr],$$

where

$$g(x)= e^{-x/2}\sum_{m=0}^{\infty} c_m m!\sum_{k=0}^{m}\frac{(-x)^k}{(k!)^2 (m-k)!}.$$

---

## Inputs and Outputs

* The partial sums in equation (2) for various $\beta$ values are read from  
  **`FIRST.txt`**, **`SECOND.txt`**, **`THIRD.txt`**, **`FOURTH.txt`**.
* The third term in equation (1) is read from **`FIFTH.txt`**.
* The final computed ground-state energy $E_0(\beta)$ is written to **`disectres.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `result.cpp`      | Main C++ source computing $E_0(\beta)$.                                |
| `FIRST.txt` … `FIFTH.txt` | Input files containing the different term expansions.            |
| `disectres.txt`   | Output file with computed ground-state energies.                         |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + GMP).            |
| `compile.job`     | SLURM job script to compile the program on an HPC cluster.              |
| `together.job`    | SLURM job script to run the program on an HPC cluster.                  |
| `run.sh`         | Legacy shell script to build and run locally on Ubuntu 22.04.           |

---

## Building and Running

Run the executable `run.sh` to compile and execute the program.  
