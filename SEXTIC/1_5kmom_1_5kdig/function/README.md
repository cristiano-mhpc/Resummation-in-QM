# Fifth Term Computation

## What this code is about

The C++ program **`function.cpp`** computes the **third term on the right-hand side** of equation (3.50):

$$E^{(3)}(\beta) = 1 + \sum_{k=0}^{\infty} (-1)^k \beta^k \mu^{-(2k+2)} + \Delta(\beta),$$

where the correction term $\Delta(\beta)$ is

$$\Delta(\beta) =\frac{\pi \beta^{(1+\nu)/2}}{\sin (\pi\nu)}\Bigl[\cos\left(\frac{\pi\nu}{2}\right)\mathrm{Im} g\left(\frac{i}{\sqrt{\beta}}\right)+ \sin\left(\frac{\pi\nu}{2}\right)\mathrm{Re} g\left(\frac{i}{\sqrt{\beta}}\right)\Bigr],$$ 

with $\nu = -\tfrac{1}{2}$ and

$$
g(x) =
e^{-x/2}\sum_{m=0}^{d} c_m m!
\sum_{k=0}^{m}
\frac{(-x)^k}{(k!)^2 (m - k)!}.
$$ 

The program reads the required $d + 1$ coefficients $c_m$ from  **`Constant.txt`** and outputs the computed values of $E^{(3)}(\beta)$   for $\beta \in \{10^{-2}, \dots, 10^{25}, 0.2, 4\}$ to **`FIFTH.txt`**.

---

## Repository Layout

| File / Folder     | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `function.cpp`    | Main C++ source computing the third term in equation (3.50).            |
| `Constant.txt`    | Input data: the d + 1 coefficients \(c_m\).                             |
| `FIFTH.txt`       | Output file containing computed values of \(E^{(3)}(\beta)\).           |
| `CMakeLists.txt`  | Modern build configuration using CMake (MPI + Boost + MPFR + MPC).      |
| `compile.job`     | SLURM script to build the project on an HPC cluster.                    |
| `together.job`    | SLURM script to run the compiled executable on an HPC cluster.          |

---

## Building and Running

To build and run the program execute the shell script `run.sh`. 
