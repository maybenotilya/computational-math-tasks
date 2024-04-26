from svd.AbstractSVD import AbstractSVD
import numpy as np
from typing import Tuple


class SimpleSVD(AbstractSVD):
    """
    Simple power iterations method
    http://www.cs.yale.edu/homes/el327/datamining2013aFiles/07_singular_value_decomposition.pdf
    """

    def __call__(self, A: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = A.astype(np.float64)
        rows, cols = A.shape

        delta = 0.1
        eps = 0.97
        lam = 2
        s = int(
            np.log(4 * np.log(2 * cols / delta) / (eps * delta)) / (2 * lam)
        )

        U = np.zeros((rows, 1))
        S = np.ndarray(0)
        V = np.zeros((1, cols))

        for _ in range(min(A.shape)):
            u, sigma, v = self._compute_single_svd(A, s)
            U = np.hstack((U, u))
            S = np.hstack((S, sigma))
            V = np.vstack((V, v))
            A -= u.dot(v) * sigma

        return U[:, 1:], S, V[1:, :]

    @staticmethod
    def _compute_single_svd(A, s) -> Tuple[np.ndarray, int, np.ndarray]:
        rows, cols = A.shape
        x = np.random.normal(0, 1, size=cols)
        B = A.T @ A

        for _ in range(s):
            x = B.dot(x)

        v = x / np.linalg.norm(x)
        sigma = np.linalg.norm(A.dot(v))
        u = A.dot(v) / sigma
        return np.reshape(u, (rows, 1)), sigma, v.reshape(1, cols)
