import numpy as np


def exact_solution(x):
    """
    Analytic solution:

        u(x) = x sin(3 pi x)
    """
    return x * np.sin(3 * np.pi * x)


def right_hand_side(x):
    """
    Right-hand side for the Poisson problem:

        -u'' = f

    For u(x) = x sin(3 pi x),

        f(x) = -6 pi cos(3 pi x) + 9 pi^2 x sin(3 pi x).
    """
    return -6 * np.pi * np.cos(3 * np.pi * x) + 9 * np.pi**2 * x * np.sin(3 * np.pi * x)


def build_matrix(N):
    """
    Build the finite difference matrix for

        -u'' = f

    with homogeneous Dirichlet boundary conditions:

        u(0) = 0, u(1) = 0.

    N is the number of subintervals, so h = 1/N.
    The unknowns are the interior values u_1, ..., u_{N-1}.
    """
    h = 1 / N
    size = N - 1

    main_diagonal = 2 * np.ones(size)
    off_diagonal = -1 * np.ones(size - 1)

    A = (
        np.diag(main_diagonal)
        + np.diag(off_diagonal, 1)
        + np.diag(off_diagonal, -1)
    )

    return A / h**2


def build_grid_and_rhs(N):
    """
    Build the interior grid points and the right-hand side vector.
    """
    h = 1 / N
    x = np.linspace(h, 1 - h, N - 1)
    b = right_hand_side(x)

    return x, b


def discrete_l2_error(x, numerical_solution):
  
    exact_values = exact_solution(x)
    error = np.sqrt(np.mean((exact_values - numerical_solution) ** 2))

    return error


def residual_norm(A, u, b):
    """
    Compute the residual norm ||Au - b||_2.
    """
    return np.linalg.norm(A @ u - b)


def direct_solver(A, b):
    """
    Solve the linear system directly using NumPy.
    """
    return np.linalg.solve(A, b)


def damped_jacobi(A, b, omega, tol=1e-10, max_iter=100000):
    """
    Damped Jacobi iteration:

        u^{k+1} = u^k + omega * D^{-1}(b - A u^k).

    Returns:
        solution, number of iterations, convergence flag.
    """
    n = len(b)
    u = np.zeros(n)

    diagonal = np.diag(A)

    for iteration in range(1, max_iter + 1):
        residual = b - A @ u
        u = u + omega * residual / diagonal

        if residual_norm(A, u, b) < tol:
            return u, iteration, True

    return u, max_iter, False


def sor(A, b, omega, tol=1e-10, max_iter=10000):
    """
    Successive Over-Relaxation method.

    The update is performed component by component.
    """
    n = len(b)
    u = np.zeros(n)

    for iteration in range(1, max_iter + 1):
        for i in range(n):
            left_sum = A[i, :i] @ u[:i]
            right_sum = A[i, i + 1:] @ u[i + 1:]

            gauss_seidel_value = (b[i] - left_sum - right_sum) / A[i, i]
            u[i] = (1 - omega) * u[i] + omega * gauss_seidel_value

        if residual_norm(A, u, b) < tol:
            return u, iteration, True

    return u, max_iter, False


def run_experiment():
    """
    Compare direct solver, damped Jacobi, and SOR
    for several grid sizes.
    """
    N_values = [8, 16, 32, 64]

    jacobi_omegas = [0.5, 0.8, 1.0]
    sor_omegas = [1.0, 1.2, 1.5, 1.8]

    print("Direct Solver, Damped Jacobi, and SOR")
    print("=====================================")
    print()

    for N in N_values:
        A = build_matrix(N)
        x, b = build_grid_and_rhs(N)

        condition_number = np.linalg.cond(A)

        print(f"N = {N}")
        print(f"h = {1 / N:.6f}")
        print(f"condition number = {condition_number:.4e}")
        print()

        # Direct solver
        direct_solution = direct_solver(A, b)
        direct_error = discrete_l2_error(x, direct_solution)

        print("Direct solver:")
        print(f"  error = {direct_error:.6e}")
        print()

        # Damped Jacobi
        print("Damped Jacobi:")
        for omega in jacobi_omegas:
            solution, iterations, converged = damped_jacobi(
                A=A,
                b=b,
                omega=omega
            )

            error = discrete_l2_error(x, solution)

            print(
                f"  omega = {omega:.1f}, "
                f"iterations = {iterations:6d}, "
                f"converged = {converged}, "
                f"error = {error:.6e}"
            )

        print()

        # SOR
        print("SOR:")
        best_omega = None
        best_iterations = None
        best_error = None

        for omega in sor_omegas:
            solution, iterations, converged = sor(
                A=A,
                b=b,
                omega=omega
            )

            error = discrete_l2_error(x, solution)

            print(
                f"  omega = {omega:.1f}, "
                f"iterations = {iterations:6d}, "
                f"converged = {converged}, "
                f"error = {error:.6e}"
            )

            if converged:
                if best_iterations is None or iterations < best_iterations:
                    best_iterations = iterations
                    best_omega = omega
                    best_error = error

        print()

        if best_omega is not None:
            print("Best SOR parameter among tested values:")
            print(f"  omega = {best_omega:.1f}")
            print(f"  iterations = {best_iterations}")
            print(f"  error = {best_error:.6e}")
        else:
            print("No tested SOR parameter converged.")

        print()
        print("-" * 60)
        print()


if __name__ == "__main__":
    run_experiment()
