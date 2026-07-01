# numerical-analysis-ii😁
This repository contains Python implementations and experiments from my Numerical Analysis II course.  Topics include numerical methods for ordinary differential equations, stability analysis, Runge-Kutta methods, implicit methods, Krylov subspace methods.

## 1. Numerical ODEs

### 1.1 ODE Stability
(1)File:

`euler_trapezoidal.py`

Comparison of Explicit Euler, Implicit Euler, and the Trapezoidal Rule for the model problem

$$
y' = \lambda y, \quad y(0)=1.
$$

This example is used to study accuracy, convergence, and stability behavior.

(2)File:

stiff_ivp_a_vs_l_stability.py

This experiment studies the stiff initial value problem

$$
y'(x)=2000(\cos x-y(x)), \quad y(0)=0, \quad x\in[0,1.5].
$$

The goal is to compare the behavior of the Implicit Euler method and the Trapezoidal Rule on a stiff problem.

### 1.2. Runge-Kutta Methods

Implementation and stability analysis of classical Runge-Kutta methods.

### 1.3. Implicit Methods

Experiments with implicit time-stepping methods such as the trapezoidal rule and Radau-type methods.

## 2. Finite Difference Methods for Boundary Value Problems

This section studies finite difference discretizations for one-dimensional boundary value problems.

### 2.1 Poisson Boundary Value Problem

We consider the Poisson boundary value problem

$$
-u'' = f \quad \text{in } (0,1),
$$

with boundary conditions

$$
u(0)=a, \qquad u(1)=b.
$$

In this experiment, we use homogeneous boundary conditions

$$
u(0)=0, \qquad u(1)=0.
$$

The right-hand side is chosen as

$$
f(x) = -6\pi \cos(3\pi x) + 9\pi^2 x \sin(3\pi x),
$$

so that the analytic solution is

$$
u(x)=x\sin(3\pi x).
$$

The finite difference approximation is used to compute numerical solutions on different grids.

File:

`bvp_convergence.py`

Experiments include:

- construction of the finite difference matrix
- comparison with the analytic solution
- computation of the discrete error
- estimation of the convergence order

### 2.2 Convergence order

The error is measured using the discrete norm

$$
|u-u_h|2 =
\left(
\frac{1}{N-1}
\sum{i=1}^{N-1}
\left(u(x_i)-u_i\right)^2
\right)^{1/2}.
$$

The convergence order is estimated using the ansatz

$$
|u-u_h|_2 \approx c h^\alpha.
$$

Here, (h) is the mesh width, (c) is a constant, and $\alpha$ is the estimated convergence order.

In the numerical experiment, the value of $\alpha$ is computed using the errors on the two finest grids.

## 3. Iterative Methods for Linear Systems

This section studies iterative methods for solving the linear systems obtained from finite difference discretizations.

The finite difference method leads to a linear system

$$
Au=b.
$$

Different solvers are used to solve the same system, and their performance is compared.

### 3.1 Direct Solver, Damped Jacobi, and SOR

File:

`jacobi_sor.py`

This experiment compares three methods:

* direct solver
* damped Jacobi method
* SOR method

For each mesh size, the script computes:

* the condition number of the finite difference matrix
* the discrete error compared with the analytic solution
* the number of iterations required by the iterative methods
* the best SOR relaxation parameter among the tested values

The damped Jacobi method is tested with

$$
\omega \in {0.5, 0.8, 1.0}.
$$

The SOR method is tested with

$$
\omega \in {1.0, 1.2, 1.5, 1.8}.
$$

The stopping criterion for the iterative methods is

$$
|Au-b|_2 < 10^{-10}.
$$

### 3.2 Damped Jacobi Method

The damped Jacobi method updates the approximation by

$$
u^{(k+1)} = u^{(k)} + \omega D^{-1}\left(b - Au^{(k)}\right).
$$

where (D) is the diagonal part of (A), and (\omega) is the damping parameter.

The goal is to observe how the damping parameter affects the number of iterations.

### 3.3 SOR Method

The SOR method, or Successive Over-Relaxation method, is based on Gauss-Seidel iteration with a relaxation parameter.

For each component, the update has the form

$$
u_i^{(k+1)} = (1-\omega)u_i^{(k)} + \omega \frac{b_i - \sum_{j<i} A_{ij}u_j^{(k+1)} - \sum_{j>i} A_{ij}u_j^{(k)}}{A_{ii}}.
$$

Here, $\omega$ is the relaxation parameter. When $\omega=1$, the method reduces to the Gauss-Seidel method.

The relaxation parameter strongly influences the convergence speed. In this experiment, several values of `omega` are tested, and the best one is chosen according to the smallest number of iterations.

### 3.4 Optimal Relaxation Parameter for SOR

File:

`optimal_sor.py`

This experiment studies the optimal relaxation parameter for the SOR method.

For the considered one-dimensional Poisson problem, the optimal relaxation parameter is given by

$$
\omega_{\mathrm{opt}}=\frac{2}{1+\sin(\pi h)}.
$$

The script compares SOR with several fixed relaxation parameters and with the theoretical optimal parameter.

The tested fixed parameters are

$$
\omega \in \{1.0, 1.2, 1.5, 1.8\}.
$$

For each mesh width, the script computes:

- the theoretical optimal relaxation parameter
- the number of SOR iterations for fixed parameters
- the number of SOR iterations using the optimal parameter
- the discrete error compared with the analytic solution

As the mesh width $h$ decreases, the value of $\omega_{\mathrm{opt}}$ approaches $2$. In the experiments, the optimal parameter usually reduces the number of iterations compared with the fixed tested values.

### 3.5 Richardson Iteration

File:

`richardson.py`

This experiment applies Richardson iteration to the linear systems obtained from the finite difference discretization of the Poisson boundary value problem.

The Richardson iteration is given by

$$
u^{(k+1)} = u^{(k)} + \omega (b - Au^{(k)}).
$$

For convergence, the damping parameter should satisfy

$$
0 < \omega < \frac{2}{\lambda_{\max}(A)}.
$$

In this experiment, a safety factor is used:

$$
\omega = 0.9 \cdot \frac{2}{\lambda_{\max}(A)}.
$$

For each mesh width, the script computes:

- the largest eigenvalue of the finite difference matrix
- the damping parameter \(\omega\)
- the number of Richardson iterations
- the discrete error compared with the analytic solution

As the mesh width \(h\) decreases, the largest eigenvalue of the finite difference matrix increases. Therefore, the damping parameter becomes smaller, and the number of iterations usually increases.

## Krylov Subspace Methods

This section contains implementations of numerical linear algebra methods based on Krylov subspaces.

A Krylov subspace generated by a matrix $A$ and a vector $B$ is

$$
\mathcal{K}_m(A,b)=\text{span}\{b, Ab, A^2b, \dots, A^{m-1}b\}.
$$

These methods are important for solving large linear systems and eigenvalue problems.

## Purpose

The goal of this repository is to organize my learning notes and coding experiments for Numerical Analysis II.😁
