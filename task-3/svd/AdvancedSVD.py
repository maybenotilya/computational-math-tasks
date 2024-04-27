from svd.AbstractSVD import AbstractSVD
import numpy as np
import scipy
from typing import Tuple
from constants import SVD_TOLERANCE


class AdvancedSVD(AbstractSVD):
    """
    Block SVD power algorithm
    https://www.emis.de/journals/ASUO/mathematics_/anale2015vol2/Bentbib_A.H.__Kanber_A..pdf
    """

    def __call__(self, A: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        rows, cols = A.shape
        k = min(rows, cols)

        V = np.random.normal(0, 1, (cols, rows))
        tol = SVD_TOLERANCE
        err = tol + 1

        while err > tol:
            Q, R = scipy.linalg.qr(A @ V, mode='economic')
            U = Q[:, :k]
            Q, R = scipy.linalg.qr(A.T @ U, mode='economic')
            V = Q[:, :k]
            S = R[:k, :k]
            err = np.linalg.norm(A @ V - U @ S)

        return U, np.diag(S), V.T
