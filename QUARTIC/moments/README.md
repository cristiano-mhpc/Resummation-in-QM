# Rayleigh–Schrödinger Perturbation Coefficients for the Quartic Anharmonic Oscillator

## What this code is about

The C++ program **`moments.cpp`** computes the Rayleigh–Schrödinger perturbation expansion coefficients $b_k$ for the ground-state energy of the quartic anharmonic oscillator:

$$E^{(2)}(\beta) = 1 + \sum_{k=1}^{\infty} b_k \beta^k.$$

The coefficients $b_k$ are computed recursively for $k = d, d-1, \dots, 1$ using:

$$b_k = \frac{3}{\sqrt{2}}\;a^{(k-1)}_{0,2}+ \frac{\sqrt{6}}{2}\;a^{(k-1)}_{0,4}, \quad k \ge 2,$$

with initial conditions

$$b_0 = 1, \qquad b_1 = \frac{3}{4}.$$

The recursion for the auxiliary coefficients is

$$a^{(k-1)}_{0,2} =-\frac{1}{4}\left[\langle 2|x^4|2\rangle a^{(k-2)}_{0,2}+ \langle 2|x^4|4\rangle a^{(k-2)}_{0,4}+ \langle 2|x^4|6\rangle a^{(k-2)}_{0,6}- \sum_{s=1}^{k-2} a^{(s)}_{0,2} b^{(k-1-s)}\right],$$

with

$$a^{(1)}_{0,2} = -\frac{\langle 2|x^4|0\rangle}{4},$$

and

$$a^{(k-1)}_{0,4} =-\frac{1}{8}\left[\langle 4|x^4|2\rangle a^{(k-2)}_{0,2}+ \langle 4|x^4|4\rangle a^{(k-2)}_{0,4}+ \langle 4|x^4|6\rangle a^{(k-2)}_{0,6}+ \langle 4|x^4|8\rangle a^{(k-2)}_{0,8}- \sum_{s=1}^{k-2} a^{(s)}_{0,4} b^{(k-1-s)}\right],$$

with

$$a^{(1)}_{0,4} = -\frac{\langle 4|x^4|0\rangle}{8}.$$

In general,

$$a^{(1)}_{0,m} = -\frac{\langle m|x^4|0\rangle}{2m}$$

and for $r \ge 2, m \ne 0$,

$$a^{(r)}_{0,m} =-\frac{1}{2m}\left[\sum_{k=1}^{\infty} \langle m|x^4|k\rangle a^{(r-1)}_{0,k} - \sum_{s=1}^{r-1} a^{(s)}_{0,m} b^{(r-s)}\right].$$

The infinite sum is finite because the only non-zero matrix elements are:

$$\begin{aligned}
\langle n|x^4|n+4\rangle &= \frac{1}{4}\sqrt{(n+1)(n+2)(n+3)(n+4)}, \\
\langle n|x^4|n+2\rangle &= \frac{1}{2}(2n+3)\sqrt{(n+1)(n+2)}, \\
\langle n|x^4|n\rangle   &= \frac{3}{2}\left(n^2 + n + \frac{1}{2}\right), \\
\langle n|x^4|n-2\rangle &= \frac{1}{2}(2n-1)\sqrt{n(n-1)}, \\
\langle n|x^4|n-4\rangle &= \frac{1}{4}\sqrt{n(n-1)(n-2)(n-3)}.
\end{aligned}
$$

All computed coefficients $b_k$ are written to **`moments.txt`**.




---

## Repository Layout

| File / Folder     | Purpose                                                                  |
|-------------------|---------------------------------------------------------------------------|
| `moments.cpp`     | C++ source computing the Rayleigh–Schrödinger perturbation coefficients. |
| `moments.txt`     | Output file containing the computed $b_k$ coefficients.                  |
| `CMakeLists.txt`  | Modern CMake build configuration (MPI + Boost + MPFR + MPC + GMP).      |
| `compile.job`     | SLURM script to compile the program on an HPC cluster.                  |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |
| `run.sh`         | To compile and run.       |

---

## Building and Running
 To build and run the program, execute the script `run.sh`. 
