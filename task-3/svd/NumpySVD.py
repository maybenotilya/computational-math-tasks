from svd.AbstractSVD import AbstractSVD
import numpy as np
from typing import Tuple


class NumpySVD(AbstractSVD):
    def __call__(self, matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        U, S, V = np.linalg.svd(matrix, full_matrices=False)
        return U, S, V
