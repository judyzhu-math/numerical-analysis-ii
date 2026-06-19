# Numerical ODE Methods

This folder contains Python implementations and numerical experiments for ordinary differential equations from my Numerical Analysis II course.

The main topics include one-step methods, error behavior, convergence, and stability for stiff and non-stiff initial value problems.

## Content
### 1. Euler and Trapezoidal Rule
File:

euler_trapezoidal.py

This experiment compares the Explicit Euler method, Implicit Euler method, and Trapezoidal Rule for the model problem

$$
y' = \lambda y, \quad y(0)=1, \quad x \in [0,1].
$$

In this example, we use

$$
\lambda = -10.
$$

The exact solution is

$$
y(x)=e^{\lambda x}.
$$

The program computes the numerical solution at (x=1) for different step numbers (N), compares it with the exact solution, and plots the absolute error.

### 2. A-stability vs. L-stability for a Stiff IVP
File:

stiff_ivp_a_vs_l_stability.py

This experiment studies the stiff initial value problem

$$
y'(x)=2000(\cos x-y(x)), \quad y(0)=0, \quad x\in[0,1.5].
$$

The goal is to compare the behavior of the Implicit Euler method and the Trapezoidal Rule on a stiff problem.

The Implicit Euler method is L-stable, so it strongly damps rapidly decaying stiff components.

The Trapezoidal Rule is A-stable but not L-stable. Therefore, it can remain stable while still producing artificial oscillations around the smooth exact solution.

This example shows why A-stability alone is not always enough for stiff initial value problems.

## Numerical Formulas

### Explicit Euler Method

$$
y_{n+1} = (1+h\lambda)y_n
$$

This method is simple, but it may become unstable when the step size is too large.

### Implicit Euler Method

$$
y_{n+1} = \frac{y_n}{1-h\lambda}
$$

This method is more stable for rapidly decaying problems.

### Trapezoidal Rule

$$
y_{n+1}=\frac{1+\frac{h\lambda}{2}}{1-\frac{h\lambda}{2}}y_n
$$

The trapezoidal rule usually has better accuracy than the Euler methods.

## Error Plot

The plot below shows how the absolute error changes as the number of steps \(N\) increases.

![Error Plot](error_plot.png)

## Explanation of the Error Plot

The plot compares the absolute error of three numerical methods:

- Explicit Euler Method
- Implicit Euler Method
- Trapezoidal Rule

The x-axis represents the number of steps \(N\). A larger \(N\) means a smaller step size

$$
h = \frac{1}{N}.
$$

The y-axis represents the absolute error at \(x=1\):

$$
\text{error} = |y_{\text{numerical}}(1) - y_{\text{exact}}(1)|.
$$

The y-axis is shown on a logarithmic scale, so a downward trend means that the error is decreasing quickly.

From the plot, we can see that as \(N\) increases, the errors of all three methods generally decrease. This means that using a smaller step size gives a more accurate numerical approximation.

The Trapezoidal Rule has the smallest error for large \(N\), showing better accuracy than the two Euler methods. The Explicit Euler Method has large errors when \(N\) is small, which shows that it can behave poorly when the step size is too large.

This experiment demonstrates the importance of both accuracy and stability in numerical methods for ordinary differential equations.



# A-stable and L-stable for a stiff IVP
## Analytic Solution

The differential equation is

$$
y'(x) = 2000(\cos x - y(x)).
$$

This can be rewritten as

$$
y'(x) + 2000y(x) = 2000\cos x.
$$

Solving this linear ODE gives

$$
y(x)=\frac{4000000}{4000001}\cos x
+\frac{2000}{4000001}\sin x
-\frac{4000000}{4000001}e^{-2000x}.
$$

The term \(e^{-2000x}\) decays extremely quickly. This fast-decaying term makes the problem stiff.

## Implicit Euler Scheme

For the implicit Euler method,

$$
y_{n+1}=y_n+h f(x_{n+1},y_{n+1}).
$$

Using

$$
f(x,y)=2000(\cos x-y),
$$

we obtain

$$
y_{n+1}=
\frac{y_n+2000h\cos(x_{n+1})}{1+2000h}.
$$

This method is L-stable, so the rapidly decaying stiff component is strongly damped.

## Trapezoidal Rule Scheme

For the trapezoidal rule,

$$
y_{n+1}=
y_n+\frac{h}{2}
\left[
f(x_n,y_n)+f(x_{n+1},y_{n+1})
\right].
$$

For this problem, this becomes

$$
y_{n+1}=
\frac{(1-1000h)y_n+1000h(\cos x_n+\cos x_{n+1})}
{1+1000h}.
$$

The trapezoidal rule is A-stable but not L-stable. Therefore, for stiff problems, it may produce artificial oscillations that are not quickly damped.

## Discussion

This example illustrates the difference between A-stability and L-stability.

Both the implicit Euler method and the trapezoidal rule are A-stable, so they can remain stable for stiff problems. However, their behavior is different.

The implicit Euler method is L-stable. This means that very stiff components are strongly damped. As a result, the numerical solution quickly follows the smooth part of the exact solution.

The trapezoidal rule is A-stable but not L-stable. Although it remains stable, it does not strongly damp the stiff transient component. As a result, artificial oscillations may appear around the smooth exact solution.

This shows that A-stability alone is not always enough for stiff initial value problems. For very stiff problems, L-stability is often preferred.
