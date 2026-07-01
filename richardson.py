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


def residual_norm(A, u, b):
    """
    Compute the residual norm ||Au - b||_2.
    """
    return np.linalg.norm(A @ u - b)


def discrete_l2_error(x, numerical_solution):
    """
    Compute the discrete L2 error compared with the analytic solution.
    """
    exact_values = exact_solution(x)
    error = np.sqrt(np.mean((exact_values - numerical_solution) ** 2))

    return error


def choose_richardson_parameter(A):
    """
    Choose a damping parameter for Richardson iteration.

    For convergence of Richardson iteration on a symmetric positive definite matrix,
    we need

        0 < omega < 2 / lambda_max(A).

    We use a safety factor 0.9.
    """
    eigenvalues = np.linalg.eigvalsh(A)
    lambda_max = np.max(eigenvalues)

    omega = 0.9 * 2 / lambda_max

    return omega, lambda_max


def richardson_iteration(A, b, omega, tol=1e-10, max_iter=100000):
    """
    Richardson iteration:

        u^{k+1} = u^k + omega * (b - A u^k).

    Returns:
        solution, number of iterations, convergence flag.
    """
    n = len(b)
    u = np.zeros(n)

    for iteration in range(1, max_iter + 1):
        residual = b - A @ u
        u = u + omega * residual

        if residual_norm(A, u, b) < tol:
            return u, iteration, True

    return u, max_iter, False


def run_experiment():
    """
    Run Richardson iteration for several mesh sizes.
    """
    N_values = [8, 16, 32, 64, 128]

    print("Richardson Iteration")
    print("====================")
    print()

    print(
        f"{'N':>8} "
        f"{'h':>12} "
        f"{'lambda_max':>16} "
        f"{'omega':>16} "
        f"{'iterations':>14} "
        f"{'converged':>12} "
        f"{'error':>16}"
    )
    print("-" * 100)

    for N in N_values:
        h = 1 / N

        A = build_matrix(N)
        x, b = build_grid_and_rhs(N)

        omega, lambda_max = choose_richardson_parameter(A)

        solution, iterations, converged = richardson_iteration(
            A=A,
            b=b,
            omega=omega,
            tol=1e-10,
            max_iter=100000
        )

        error = discrete_l2_error(x, solution)

        print(
            f"{N:8d} "
            f"{h:12.6f} "
            f"{lambda_max:16.6e} "
            f"{omega:16.6e} "
            f"{iterations:14d} "
            f"{str(converged):>12} "
            f"{error:16.6e}"
        )


if __name__ == "__main__":
    run_experiment()
