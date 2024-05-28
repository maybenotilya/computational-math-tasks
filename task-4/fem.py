import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import quad

"""
This solver only works for Sturm–Liouville problem below: 
    y'' - λy = -2λsin(sqrt(λ)x)
    y(0) = y(n * pi) = 0

Sturm–Liouville problem looks like:
    (p(x)y')' - q(x)y = -f(x), 0 <= x <= l
    y(0) = y(l) = 0

For this case:
    p(x) = 1
    q(x) = λ
    f(x) = 2λsin(sqrt(λ)x)
Real solution is:
    y = sin(sqrt(λ)x)
For this solver we assume that x[j] - x[j - 1] = h for all j = 1..N.

As sin(x) = 0 iff x = n * pi, n in Z, sqrt(λ)n must be integer. For example, it will always be integer if λ = m^2, m in Z.

"""

def phi(j, x, x_knots):
    """
    Phi represents family of basis functions for FEM solver.
    Args:
        x: Argument of function
        j: Index of basis fuction
        x_knots: Array of x values in knots of the grid
    """
    N = len(x_knots) - 1
    h = x_knots[1] - x_knots[0]
    match j:
        case 0:
            return (x_knots[1] - x) / h if x_knots[0] <= x <= x_knots[1] else 0
        case k if k == N:
            return (x - x_knots[N - 1]) / h if x_knots[N - 1] <= x <= x_knots[N] else 0
        case _:
            if x_knots[j - 1] <= x <= x_knots[j]:
                return (x - x_knots[j - 1]) / h
            if x_knots[j] <= x <= x_knots[j + 1]:
                return (x_knots[j + 1] - x) / h
            return 0

def u(x, x_knots, y_knots):
    """
      Function u(x) represents solution of equation.

      Args:
        x: Value of x 
        x_knots: Array of x values in knots of the grid.
        y_knots: Array of y values in knots of the grid.

      Returns:
        Value of approximate solution in point x.
    """
    l, r = 0, len(y_knots) - 1
    while r - l > 1:
        m = (l + r) // 2
        if x > x_knots[m]:
            l = m
        else:
            r = m
    return y_knots[l] * phi(l, x, x_knots) + y_knots[r] * phi(r, x, x_knots)


def thomas_solver(a, b, c, d):
    """
      Solve tridioganal system of linear equations using Thomas method.
      https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm

      Args:
        a: Elements below main diagonal.
        b: Elements on main diagonal.
        c: Elements above main diagonal.
        d: Elements on the right side of equation.

      Returns:
        Solution of system
    """
    n = len(b)
    y = np.zeros(n)

    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]

    y[n - 1] = d[n - 1] / b[n - 1]

    for i in range(n - 2, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    return y

def dot_phi_A(x_knots, λ, j, k):
    """
      Finds dot production (phi_j, phi_k)_A in D_A.
      D_A = { u in C^2[0, l] : u(0) = u(l) = 0 }
      (u, v)_A = Int_0^l [p(x)u'(x)v'(x) + q(x)u(x)v(x)]dx,   u, v in D_A

      Args:
        x_knots: Array of x values in knots of the grid.
        λ: Lambda value for equation above.
        j, k: Indeces of basis functions phi.

      Returns:
        Scalar value of dot product (phi_j and phi_k)_A.
    """
    h = x_knots[1] - x_knots[0]
    # As dot product is symmetric
    j, k = min(j, k), max(j, k)
    if j == k:
        return (
            x_knots[j + 1] 
            + λ * x_knots[j] ** 2 * x_knots[j + 1] 
            - λ * x_knots[j] * x_knots[j + 1] ** 2 
            + (λ * x_knots[j + 1] ** 3) / 3 
            - x_knots[j - 1] 
            - λ * x_knots[j] ** 2 * x_knots[j - 1] 
            + λ * x_knots[j] * x_knots[j - 1] ** 2 
            - (λ * x_knots[j - 1] ** 3) / 3
                ) / (h ** 2)
    if j + 1 == k:
        return (6 - λ * (x_knots[j] - x_knots[j - 1]) ** 2) * (x_knots[j] - x_knots[j + 1]) / (6 * h ** 2)

def dot_f_phi(x_knots, λ, k):
    """
      Finds dot production (f, phi_k).
      (u, v) = Int_0^l u(x)v(x)dx,   u, v in L_2(0, l)

      Args:
        x_knots: Array of x values in knots of the grid.
        λ: Lambda value for equation above.
        k: Index of basis functions phi.

      Returns:
        Scalar value of dot product (f, phi_k).
    """
    h = x_knots[1] - x_knots[0]

    return 2 * (
        - np.sqrt(λ) * (x_knots[k] - x_knots[k + 1]) * np.cos(np.sqrt(λ) * x_knots[k]) 
        + np.sin(np.sqrt(λ) * x_knots[k]) 
        - np.sin(np.sqrt(λ) * x_knots[k + 1])
        - np.sqrt(λ) * (x_knots[k] - x_knots[k - 1]) * np.cos(np.sqrt(λ) * x_knots[k]) 
        + np.sin(np.sqrt(λ) * x_knots[k]) 
        - np.sin(np.sqrt(λ) * x_knots[k - 1])
    ) / h

def build_equation(x_knots, λ, N):
    """
      Build tridioganal matrix for equation above.

      Args:
        x_knots: Array of x values in knots of the grid.
        λ: Lambda value for equation above.
        N: Number of knots

      Returns:
        Tuple (a, b, c, d)
        a: Elements below main diagonal.
        b: Elements on main diagonal.
        c: Elements above main diagonal.
        d: Elements on the right side of equation.
    """
    a = np.zeros(N - 1)
    b = np.zeros(N - 1)
    c = np.zeros(N - 1)
    d = np.zeros(N - 1)

    for i in range(1, N):
        if i >= 2:
            a[i - 1] = dot_phi_A(x_knots, λ, i - 1, i)
        if i + 1 < N:
            c[i - 1] = dot_phi_A(x_knots, λ, i + 1, i)
        b[i - 1] = dot_phi_A(x_knots, λ, i, i)
        d[i - 1] = dot_f_phi(x_knots, λ, i)

    return a, b, c, d

def get_y_knots(x_knots, λ, N):
    a, b, c, d = build_equation(x_knots, λ, N)
    y = thomas_solver(a, b, c, d)
    # Boundary conditions
    y = np.hstack((np.zeros(1), y, np.zeros(1)))
    return y

def solution(x, λ):
    """
    Represents real solution y = sin(sqrt(λ)x)
    """
    return np.sin(np.sqrt(λ) * x)

def plot_solution(x_knots, y_knots, λ, name='plot.png'):
    """
    Plots two solutions: real and approximate ones
    """
    l, r = 0, x_knots[-1]
    x_axis = np.linspace(l, r, 10_000)
    fig, ax = plt.subplots(figsize=(18, 7))
    ax.set_title(f'N = {len(x_knots) - 1}, λ = {λ}')
    ax.plot(x_axis, solution(x_axis, λ), label='real', linewidth=0.7)
    ax.plot(x_knots, y_knots, label='approximate', linewidth=0.7)
    ax.legend()
    fig.savefig(name)
    plt.close(fig)

def plot_errors(errors, hs, name='errors.png'):
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_title(f'Correlation between h and Error')
    ax.loglog(hs, errors, 'o')
    ax.set_xlabel('h')
    ax.set_ylabel('Error')
    ax.grid(True)
    fig.savefig(name)
    plt.close(fig)

def compute_error(x_knots, y_knots, λ):
    """
    Compute maximum value of |u(x) - y(x)|
    Args:
        x_knots: Array of x values in knots of the grid.
        y_knots: Array of y values in knots of the grid.
    Returns:
        Maximum value of error
    """
    max_error = 0
    h = (x_knots[1] - x_knots[0]) / 10
    for i in range(10 * (len(x_knots) - 1) + 1):
        err = abs(u(i * h, x_knots, y_knots) - solution(i * h, λ))
        if err > max_error:
            max_error = err
    
    return max_error

def normL2(f, l):
    """
    Calculates ||f||_{L_2(r, l)}
    """
    norm, _ = np.sqrt(quad(lambda x: f(x) ** 2, 0, l, limit=1000))
    return norm

def f(x, λ):
    """
    Represents f(x) = 2λsin(sqrt(λ)x) above
    """
    return 2 * λ * np.sin(np.sqrt(λ) * x)

def c(λ, l):
    """
    Computes c = λl^2 / 4 + 1
    """
    return λ * l ** 2 / 4 + 1

def c_stroke(λ, l):
    """
    Computes c' = 4 * pi * sqrt(1 + λl^2 / 4).
    J_M = l as we use x(e) = le, x(0) = 0, x(1) = l
    """
    return l * np.sqrt(1 + λ * l ** 2 / 4)

def estimate_error(λ, h, l, x_knots, y_knots):
    """
    Estimages error of FEM.
    Returns:
        Difference between ||y - y_h|| and (c'c)^2 * h^2 * ||f||.
        If difference is negative, that proves theoretical coverage speed of O(h^2) 
    """
    return normL2(lambda x: solution(x, λ) - u(x, x_knots, y_knots), l) - (c(λ, l) * c_stroke(λ, l)) ** 2 * h ** 2 * normL2(lambda x: f(x, λ), l)


if __name__ == '__main__':
    Ns = [10, 50, 100, 500, 1000, 5000, 10_000]
    λs = [0.25, 1, 4, 16, 25, 100]
    l =  2 * np.pi
    for N in Ns:
        for λ in λs:
            x_knots = np.linspace(0, l, N + 1)
            y_knots = get_y_knots(x_knots, λ, N)
            if λ == 4:
                plot_solution(x_knots, y_knots, λ, f'./plots/N{N}_λ{λ}.png')
            err = compute_error(x_knots, y_knots, λ)
            h = x_knots[1] - x_knots[0]
            estimate = estimate_error(λ, h, l, x_knots, y_knots)
            print(f'N = {N: >5} | λ = {λ: >4} | Estimate = {estimate:.4f}')