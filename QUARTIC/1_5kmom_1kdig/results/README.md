# Ground-State Energy of the Quartic Anharmonic Oscillator  
(MPI + Boost.MPFR/GMP)

## What this code is about

The C++ code **`result.cpp`** computes the ground–state energy  
$E_0(\beta)$ of the **quartic anharmonic oscillator**:

$$E(\beta) = -\frac{1}{2}\left[1 -\sum_{k=0}^{\infty} (-1)^k \mu^{-(k+1)} \beta^k-\frac{\pi g(-1/\beta) \beta^{2/3}}{\sin(2\pi/3)}\right].$$

---

### Second term of equation (1)

The second term is expanded as
$$\sum_{k=0}^{\infty} (-1)^k \mu^{-(k+1)} \beta^k=\sum_{k=0}^{d} (-1)^k \beta^k (A_k + B_k + C_k)+\sum_{k=d+1}^{\infty} (-1)^k \beta^k D_k.$$

with
$$A_k =\sum_{m=0}^{k} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m - l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx,$$

$$B_k =\sum_{m=k+1}^{d} c_m m!\sum_{l=0}^{k}\frac{(-1)^l}{(l!)^2 (m - l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx,$$

$$C_k =\sum_{m=k+1}^{d} c_m m!\sum_{l=k+1}^{m}\frac{(-1)^l \Gamma(l - k - \nu) 2^{l-k-\nu}}{(l!)^2 (m - l)!},$$

$$D_k =\sum_{m=0}^{d} c_m m!\sum_{l=0}^{m}\frac{(-1)^l}{(l!)^2 (m - l)!}\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx.$$

The finite-part integrals are
$$\int_{0}^{\infty} e^{-x/2} x^{k+\nu+1-l}dx=(-1)^{k-l+1}\left(\frac{1}{2}\right)^{k-l+1+\nu}\pi \Gamma(k - l + 1 + \nu)\sin(\pi \nu),$$
with $\nu = \tfrac{2}{3}$.

---

### Third term of equation (1)

The function $g(x)$ is
$$g(x) =e^{-x/2}\sum_{m=0}^{d} c_m m!\sum_{k=0}^{m}\frac{(-x)^k}{(k!)^2 (m - k)!}.$$

---

## Inputs and Outputs

* The partial sums in equation (2) for various \(\beta\) values are read from  
  **`FIRST.txt`**, **`SECOND.txt`**, **`THIRD.txt`**, **`FOURTH.txt`**.
* The third term in equation (1) is read from **`FIFTH.txt`**.
* The final computed ground-state energy \(E_0(\beta)\) is written to **`disectres.txt`**.

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
