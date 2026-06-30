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

    with homogeneous Dirichlet boundary conditions.
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
    Build the interior grid points and right-hand side vector.
    """
    h = 1 / N
    x = np.linspace(h, 1 - h, N - 1)
    b = right_hand_side(x)

    return x, b


def residual_norm(A, u, b):
    """
    Compute ||Au - b||_2.
    """
    return np.linalg.norm(A @ u - b)


def discrete_l2_error(x, numerical_solution):
    """
    Compute the discrete L2 error compared with the analytic solution.
    """
    exact_values = exact_solution(x)
    error = np.sqrt(np.mean((exact_values - numerical_solution) ** 2))

    return error


def optimal_sor_parameter(h):
    """
    Compute the optimal SOR relaxation parameter for this model problem:

        omega_opt = 2 / (1 + sin(pi h)).
    """
    return 2 / (1 + np.sin(np.pi * h))


def sor(A, b, omega, tol=1e-10, max_iter=10000):
    """
    Successive Over-Relaxation method.
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


def compare_sor_parameters():
    """
    Compare SOR with several fixed omega values and with the theoretical optimal omega.
    """
    N_values = [8, 16, 32, 64, 128]
    fixed_omegas = [1.0, 1.2, 1.5, 1.8]

    print("Optimal SOR Relaxation Parameter")
    print("================================")
    print()

    for N in N_values:
        h = 1 / N
        A = build_matrix(N)
        x, b = build_grid_and_rhs(N)

        omega_opt = optimal_sor_parameter(h)

        print(f"N = {N}")
        print(f"h = {h:.6f}")
        print(f"omega_opt = {omega_opt:.6f}")
        print()

        print("Fixed omega values:")

        best_fixed_omega = None
        best_fixed_iterations = None

        for omega in fixed_omegas:
            solution, iterations, converged = sor(A, b, omega)
            error = discrete_l2_error(x, solution)

            print(
                f"  omega = {omega:.1f}, "
                f"iterations = {iterations:6d}, "
                f"converged = {converged}, "
                f"error = {error:.6e}"
            )

            if converged:
                if best_fixed_iterations is None or iterations < best_fixed_iterations:
                    best_fixed_iterations = iterations
                    best_fixed_omega = omega

        print()

        print("Theoretical optimal omega:")

        solution, iterations, converged = sor(A, b, omega_opt)
        error = discrete_l2_error(x, solution)

        print(
            f"  omega_opt = {omega_opt:.6f}, "
            f"iterations = {iterations:6d}, "
            f"converged = {converged}, "
            f"error = {error:.6e}"
        )

        print()

        if best_fixed_omega is not None:
            print("Comparison:")
            print(f"  best tested fixed omega = {best_fixed_omega}")
            print(f"  iterations with best fixed omega = {best_fixed_iterations}")
            print(f"  iterations with omega_opt = {iterations}")

        print()
        print("-" * 70)
        print()


if __name__ == "__main__":
    compare_sor_parameters()
