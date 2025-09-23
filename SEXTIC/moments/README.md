# Moments – Rayleigh–Schrödinger Perturbation Expansion

## What this code is about

The C++ program **`moments.cpp`** computes the Rayleigh–Schrödinger perturbation expansion coefficients $b_k$  
for the ground-state energy of the sextic anharmonic oscillator:

$$
E(\beta) = 1 + \sum_{k=1}^{\infty} b_k \, \beta^k
$$

up to orders $k = d, d-1, \dots ,1$, using the recursive formulas

$$
b_0 = 1, \qquad
b_1 = \frac{15}{8}
$$

$$
b_k = \frac{54\sqrt{2}}{8} a^{(k-1)}_{0,2}
     + \frac{15\sqrt{4!}}{8} a^{(k-1)}_{0,4}
     + \frac{\sqrt{6!}}{8} a^{(k-1)}_{0,6}, \quad k \ge 2
$$ 

where the auxiliary coefficients $a^{(r)}_{0,m}$ satisfy

$$
a^{(0)}_{0,m} = 0, \qquad
a^{(1)}_{0,m} = -\frac{\langle m|x^6|0\rangle}{2m}
$$


$$
a^{(r)}_{0,m} = -\frac{1}{2m}
\left[
\sum_{k=1}^\infty \langle m|x^6|k\rangle a^{(r-1)}_{0,k} - \sum_{s=1}^{r-1} a^{(s)}_{0,m}\, b_{r-s}
\right], \quad r \ge 2, \; m \neq 0
$$

Only a finite set of matrix elements $\langle n|x^6|n+\ell\rangle$ are non-zero:

$$
\begin{aligned}
\langle n|x^6|n+6\rangle &= \frac{1}{8}\sqrt{(n+6)(n+5)(n+4)(n+3)(n+2)(n+1)} \\
\langle n|x^6|n+4\rangle &= \frac{3}{8}\sqrt{(n+4)(n+3)(n+2)(n+1)(2n+5)} \\
\langle n|x^6|n+2\rangle &= \frac{15}{8}(n^2 + 3n + 3)\sqrt{(n+1)(n+2)} \\
\langle n|x^6|n\rangle   &= \frac{5}{8}\bigl(4n^3 + 6n^2 + 8n + 3\bigr) \\
\langle n|x^6|n-2\rangle &= \frac{15}{8}\sqrt{n(n-1)(n^2 - n + 1)} \\
\langle n|x^6|n-4\rangle &= \frac{3}{8}\sqrt{n(n-1)(n-2)(n-3)(2n-3)} \\
\langle n|x^6|n-6\rangle &= \frac{1}{8}\sqrt{n(n-1)(n-2)(n-3)(n-4)(n-5)}
\end{aligned}
$$ 

The computed coefficients $b_k$ are written to the file **`moments.txt`**.

---

## Repository Layout

| File / Folder         | Purpose                                                                 |
|-----------------------|-------------------------------------------------------------------------|
| `moments.cpp`         | Main C++ source computing the coefficients \(b_k\).                     |
| `moments.txt`         | Output file containing computed coefficients.                           |
| `CMakeLists.txt`      | Build configuration for CMake (MPI + Boost + MPFR).                     |
| `run.sh`              | Portable run script: configures, builds, and runs (local or SLURM).     |
| `mpfr.sh`             | Legacy local build/run script using manual `mpicxx` calls.             |
| `compile.job`         | SLURM script to build on an HPC cluster.                                |
| `together.job`        | SLURM script to run the program on an HPC cluster.                      |

---

## Building and Running

### 1. Local Ubuntu (with CMake)

To build and run the program simple run `./run.sh` 

